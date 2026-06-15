/*!
   \file do_barrier.cpp
   \brief Creates barriers
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2026
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <cstddef>
#include <assert.h>

#include <iostream>
using std::endl;
using std::ios;

#include <array>
using std::array;

#include "cme.h"
#include "simulation.h"
#include "coast_landform.h"

//===============================================================================================================================
//! Simulates barrier formation
//===============================================================================================================================
int CSimulation::nDoBarrierFormation(void)
{
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         // If waves are off-shore, then do nothing, move to the next coast point
         if (! m_VCoast[nCoast].bGetWavesOnShore(nCoastPoint))
            continue;

         // OK, waves are on-shore
         CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint);
         int nCoastLandform = pCoastLandform->nGetLandFormCategory();

         // If this isn't a beach then do nothing, move to the next coast point
         if (nCoastLandform != LF_DRIFT_BEACH)
            continue;

         // OK it is a beach. Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

         // And get the this-iteration runup for this coast point
         double const dRunUp = m_VCoast[nCoast].dGetRunUp(nCoastPoint);

         // Calc total wave elevation
         double const dWaveElev = m_dThisIterSWL + dRunUp;

         // TODO calculate inland movement of sand and gravel
         bool bWavesDownCoast = m_VCoast[nCoast].bGetWavesDownCoast(nCoastPoint);
      }
   }

   int nRet = nMoveUnconsLandward();
   if (nRet != RTN_OK)
      return nRet;

   return RTN_OK;
}


//===============================================================================================================================
//! Uses runup to calculate landward movement of sand and gravel unconsolidarted sediment
//===============================================================================================================================
int CSimulation::nMoveUnconsLandward(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Calculating cliff collapse" << endl;

   int nRet;

   // First go along each coastline and at each point on the coastline, move sand and gravel unconsolidated sediment landward
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {/*
         CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint);

         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

         // Is there some talus protecting this cell?
         double dTalusDepth = m_pRasterGrid->m_Cell[nX][nY].dGetAllTalusDepth();
         double dInvTalusProtection = 1;
         if (dTalusDepth > 0)
         {
            // TEST TODO Assume a linear relationship, with minimum value 0.5
            double dCliffHeightAboveSWL = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() - m_dThisIterSWL;
            dInvTalusProtection = tMin(tMax((dTalusDepth / dCliffHeightAboveSWL), 1.0), 0.5);

            LogStream << m_ulIter << ":\t cell[" << nX << "][" << nY << "] talus depth = " << dTalusDepth << " cliff height (in talus) above SWL = " << dCliffHeightAboveSWL << " inverse talus protection factor = " << dInvTalusProtection << endl;
         }

         // First get wave energy for the coastal landform object
         double const dWaveHeightAtCoast = m_VCoast[nCoast].dGetCoastWaveHeight(nCoastPoint);

         // If the waves at this point are off-shore, then do nothing, just move to next coast point
         if (bFPIsEqual(dWaveHeightAtCoast, DBL_NODATA, TOLERANCE))
            continue;

         // OK we have on-shore waves so get the previously-calculated wave energy
         double const dWaveEnergy = m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nCoastPoint) * dInvTalusProtection;

         // And save the accumulated value
         pCoastLandform->IncTotAccumWaveEnergy(dWaveEnergy);

         int const nCat = pCoastLandform->nGetLandFormCategory();

         // Is this a cliff?
         if (nCat == LF_CLIFF)
         {
            // It is, so get the cliff object
            CRWCliff* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

            // And do the notch incision, if any. Note that we consider sediment eroded due to notch incision to be still in place until cliff collapse, i.e. the sediment which filled the notch, pre-incision, is assumed to remain there. If the notch is eventually incised sufficiently to cause cliff collapse, then the sediment from the notch volume is included with the above-notch talus
            if (! bIncreaseCliffNotchIncision(nCoast, nX, nY, pCliff, dWaveEnergy))
               // No incision of this cliff
               continue;

            // OK, we've had some incision of this coast cliff. So is the notch now incised enough to cause collapse?  (either because the overhang is greater than the threshold overhang, or because there is no sediment remaining)?
            if (pCliff->bReadyToCollapse(m_dNotchIncisionAtCollapse))
            {
               // It is ready to collapse, so do the cliff collapse
               nRet = nDoCliffCollapse(nCoast, nX, nY, pCliff->dGetNotchApexElev());
               if (nRet != RTN_OK)
               {
                  if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  {
                     LogStream << m_ulIter << ":\t" << WARN << "problem with coast cliff collapse, continuing however" << endl;

                     if (nRet == RTN_ERR_CLIFF_NOT_IN_POLYGON)
                        LogStream << m_ulIter << ":\t coast cliff-collapse cell not in a polygon" << endl;
                     else if (nRet == RTN_ERR_CLIFF_NOTCH)
                        LogStream << m_ulIter << ":\t coast cliff notch is incised into basement" << endl;
                     else if (nRet == RTN_ERR_NO_TOP_LAYER_DURING_CLIFF_COLLAPSE_CALC)
                        LogStream << m_ulIter << ":\t no top layer during coast cliff collapse" << endl;
                  }
               }

               pCliff->SetCliffCollapsed();
            }
         }*/
      }
   }

   return RTN_OK;
}
