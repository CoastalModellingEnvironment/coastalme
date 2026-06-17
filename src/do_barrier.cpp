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

#include <ios>
using std::fixed;
using std::scientific;

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
      int const nCoastLen = m_VCoast[nCoast].nGetCoastlineSize();

      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         // If waves are off-shore, then do nothing
         if (! m_VCoast[nCoast].bGetWavesOnShore(nCoastPoint))
            continue;

         // If we are not in the active zone, then do nothing
         double dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint);
         if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
            continue;

         // OK, waves are on-shore,and we are in the active zone
         CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint);
         int nCoastLandform = pCoastLandform->nGetLandFormCategory();

         // If this isn't a beach then do nothing
         if (nCoastLandform != LF_DRIFT_BEACH)
            continue;

         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nCoastX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         int const nCoastY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

         // Get the top layer
         int nTopLayer = m_pRasterGrid->m_Cell[nCoastX][nCoastY].nGetTopNonZeroLayerAboveBasement();

         // Safety check
         if ((nTopLayer == NO_NONZERO_THICKNESS_LAYERS) || (nTopLayer == INT_NODATA))
            continue;

         // Any uncons sand or uncons coarse here?
         double dSand = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
         double dCoarse = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();

         if ((bFPIsEqual(dSand, 0.0, TOLERANCE)) && (bFPIsEqual(dCoarse, 0.0, TOLERANCE)))
            //  No uncons sand or gravel so do nothing
            continue;

         // If there is talus on this cell then do nothing
         if (m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->bHasTalus())
            continue;

         // OK so far, so now try to move some uncons sand or uncons gravel inland
         double dCoastElev = m_pRasterGrid->m_Cell[nCoastX][nCoastY].dGetAllSedTopElevIncTalus();

         // Get the coastline points before and after this one
         int nCoastXBefore = nCoastX;
         int nCoastYBefore = nCoastY;
         int nCoastXAfter = nCoastX;
         int nCoastYAfter = nCoastY;

         if (nCoastPoint > 0)
         {
            nCoastXBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint - 1)->nGetX();
            nCoastYBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint - 1)->nGetY();
         }

         if (nCoastPoint < nCoastLen - 1)
         {
            nCoastXAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint + 1)->nGetX();
            nCoastYAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint + 1)->nGetY();
         }

         // Get the this-iteration runup for this coast point
         double const dRunUp = m_VCoast[nCoast].dGetRunUp(nCoastPoint);

         // Calc total elevation of runup
         // double const dRunUpTopElev = m_dThisIterSWL + dRunUp;
         double const dRunUpTopElev = m_dThisIterMeanSWL + dRunUp;

         // TODO calculate inland movement of sand and gravel
         int nHanded = m_VCoast[nCoast].nGetSeaHandedness();      // RH = 0, LH = 1
         int const nCoastHand = m_VCoast[nCoast].nGetSeaHandedness();

         bool bWavesDownCoast = m_VCoast[nCoast].bGetWavesDownCoast(nCoastPoint);
         double dTangentToCoast = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);
         double dBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nCoastPoint);

         // // TEST
         // double dDiff = dKeepWithin360(dBreakingWaveAngle - dTangentToCoast + 180);
         //
         // double dDummy = -1;
         // if ((dDiff >= 315) || (dDiff < 45))
         //    // Inland diagonal upcoast
         //    dDummy = 1;
         // else if ((dDiff >= 45) && (dDiff < 135))
         //    // Inland 90 degrees
         //    dDummy = 2;
         // else if ((dDiff >= 135) && (dDiff < 225))
         //    // Inland diagonal downcoast
         //    dDummy = 3;
         // else if ((dDiff >= 225) && (dDiff < 315))
         //    // Inland 90 degrees
         //    dDummy = 4;
         //
         // LogStream << "dDummy = " << dDummy << endl;

         bool bIsACellLessThanRunupElev = true;
         int n = 0;
         double dLastSandMoved = DBL_NODATA;
         double dLastCoarseMoved = DBL_NODATA;
         CGeom2DIPoint PtiLast(INT_NODATA, INT_NODATA);

         do
         {
            n++;
            CGeom2DIPoint const PtiTmp = PtiGetPerpendicular(nCoastXBefore, nCoastYBefore, nCoastXAfter, nCoastYAfter, n * m_dCellSide, nCoastHand);

            // Safety check
            if (! bIsWithinValidGrid(&PtiTmp))
               break;

            // Prevent duplication due to rounding
            if (PtiTmp == PtiLast)
            {
               PtiLast = PtiTmp;
               continue;
            }

            int nTmpX = PtiTmp.nGetX();
            int nTmpY = PtiTmp.nGetY();
            double dCellElev = m_pRasterGrid->m_Cell[nTmpX][nTmpY].dGetAllSedTopElevIncTalus();

            if (dCellElev < dRunUpTopElev)
            {
               bIsACellLessThanRunupElev = true;

               if (PtiLast.nGetX() == INT_NODATA)
               {
                  PtiLast.SetX(nCoastX);
                  PtiLast.SetY(nCoastY);
               }

               // (runuptop - elev) / (runuptop - coastelev)
               double dWeight = (dRunUpTopElev - dCellElev) / (dRunUpTopElev - dCoastElev);

               // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t possible barrier uncons movement, coast point = " << nCoastPoint << " [" << nCoastX << "][" << nCoastY << "] = {" << dGridCentroidXToExtCRSX(nCoastX) << ", " << dGridCentroidYToExtCRSY(nCoastY) << "}, last point [" << PtiLast.nGetX() << "][" << PtiLast.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiLast.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiLast.nGetY()) << "}, this point [" << nTmpX << "][" << nTmpY << "] = {" << dGridCentroidXToExtCRSX(nTmpX) << ", " << dGridCentroidYToExtCRSY(nTmpY) << "} dRunUp = " << dRunUp << " dRunUpTopElev = " << dRunUpTopElev << " cell elev = " << dCellElev << " dWeight = " << dWeight << endl;

               int nRet = nMoveUnconsLandward(&PtiLast, &PtiTmp, dWeight, dLastSandMoved, dLastCoarseMoved);
               // if (nRet != RTN_OK)
                  // return nRet;

               // Don't continue if last sediment moved was a tiny amount
               if (dLastSandMoved < SED_ELEV_TOLERANCE)
                  break;

               if (dLastCoarseMoved < SED_ELEV_TOLERANCE)
                  break;
            }
            else
               bIsACellLessThanRunupElev = false;

            PtiLast = PtiTmp;

         } while (bIsACellLessThanRunupElev);
      }
   }

   return RTN_OK;
}


//===============================================================================================================================
//! Uses runup to calculate landward movement of sand and gravel unconsolidated sediment
//===============================================================================================================================
int CSimulation::nMoveUnconsLandward(CGeom2DIPoint const* pPtiFrom, CGeom2DIPoint const* pPtiTo, double const dWeight, double& dLastSandMoved, double& dLastCoarseMoved)
{
   int nXFrom = pPtiFrom->nGetX();
   int nYFrom = pPtiFrom->nGetY();
   int nXTo = pPtiTo->nGetX();
   int nYTo = pPtiTo->nGetY();

   int nTopLayer = m_pRasterGrid->m_Cell[nXFrom][nYFrom].nGetTopNonZeroLayerAboveBasement();

   // Safety check
   if ((nTopLayer == NO_NONZERO_THICKNESS_LAYERS) || (nTopLayer == INT_NODATA))
   {
      LogStream << "Down to basement" << endl;
      return RTN_ERR_BASEMENT_DURING_BARRIER_CREATION;
   }

   double dSandThis;
   double dCoarseThis;

   if (bFPIsEqual(dLastSandMoved, DBL_NODATA, TOLERANCE))
      dSandThis = m_pRasterGrid->m_Cell[nXFrom][nYFrom].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
   else
      dSandThis = dLastSandMoved;

   if (bFPIsEqual(dLastCoarseMoved, DBL_NODATA, TOLERANCE))
      dCoarseThis = m_pRasterGrid->m_Cell[nXFrom][nYFrom].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
   else
      dCoarseThis = dLastCoarseMoved;

   // TEST
   double dFractSand = 0.8;
   double dFractCoarse = 0.4;

   double dSandToMove = 0;
   double dCoarseToMove = 0;

   if (dSandThis > 0)
   {
      dSandToMove = dSandThis * dFractSand * dWeight;

      // Don't move tiny amounts
      if (dSandToMove > SED_ELEV_TOLERANCE)
      {
         m_pRasterGrid->m_Cell[nXFrom][nYFrom].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandThis - dSandToMove);
         m_pRasterGrid->m_Cell[nXTo][nYTo].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandToMove);

         // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t  barrier sand movement, from [" << nXFrom << "][" << nYFrom << "] = {" << dGridCentroidXToExtCRSX(nXFrom) << ", " << dGridCentroidYToExtCRSY(nYFrom) << "} to [" << nXTo << "][" << nYTo << "] = {" << dGridCentroidXToExtCRSX(nXTo) << ", " << dGridCentroidYToExtCRSY(nYTo) << "} sand depth moved = " << scientific << dSandToMove << fixed << endl;
      }
   }

   if (dCoarseThis > 0)
   {
      dCoarseToMove = dCoarseThis * dFractCoarse * dWeight;

      if (dCoarseToMove > SED_ELEV_TOLERANCE)
      {
         m_pRasterGrid->m_Cell[nXFrom][nYFrom].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseThis - dCoarseToMove);
         m_pRasterGrid->m_Cell[nXTo][nYTo].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseToMove);

         // if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
            LogStream << m_ulIter << ":\t  barrier coarse movement, from [" << nXFrom << "][" << nYFrom << "] = {" << dGridCentroidXToExtCRSX(nXFrom) << ", " << dGridCentroidYToExtCRSY(nYFrom) << "} to [" << nXTo << "][" << nYTo << "] = {" << dGridCentroidXToExtCRSX(nXTo) << ", " << dGridCentroidYToExtCRSY(nYTo) << "} coarse depth moved = " << scientific << dCoarseToMove << fixed << endl;
      }
   }

   // For next time
   dLastSandMoved = dSandToMove;
   dLastCoarseMoved = dCoarseToMove;

   return RTN_OK;
}
