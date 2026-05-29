/*!
   \file update_grid.cpp
   \brief Updates the raster grid
   \details TODO 001 A more detailed description of this routine.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <iostream>
using std::endl;

#include <cfloat>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cme.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! At the end of each timestep, updates all cells in the raster grid and does some per-timestep accounting
//===============================================================================================================================
int CSimulation::nEndOfTimestepUpdateGrid(void)
{
   // Go through all cells in the raster grid and calculate some this-timestep totals
   m_dThisIterTopElevMax = -DBL_MAX;
   m_dThisIterTopElevMin = DBL_MAX;

   // Resize vectors to store per-polygon totals of fine, sand, and coarse talus
   int nCoastSize = static_cast<int>(m_VCoast.size());
   m_VdFineTalus.resize(nCoastSize);
   for (int n = 0; n < (nCoastSize); n++)
   {
      int nPoly = m_VCoast[n].nGetNumPolygons();
      m_VdFineTalus[n].resize(nPoly);
   }

   m_VdSandTalus.resize(nCoastSize);
   for (int n = 0; n < (nCoastSize); n++)
   {
      int nPoly = m_VCoast[n].nGetNumPolygons();
      m_VdSandTalus[n].resize(nPoly);
   }

   m_VdCoarseTalus.resize(nCoastSize);
   for (int n = 0; n < (nCoastSize); n++)
   {
      int nPoly = m_VCoast[n].nGetNumPolygons();
      m_VdCoarseTalus[n].resize(nPoly);
   }

   // Initialize reduction variables to zero
   m_ulThisIterNumCoastCells = 0;
   m_dThisIterTotSeaDepth = 0;

   // Use OpenMP parallel reduction for thread-safe accumulation and min/max calculations
#ifdef _OPENMP
#pragma omp parallel for collapse(2)                                 \
    reduction(+ : m_ulThisIterNumCoastCells, m_dThisIterTotSeaDepth) \
    reduction(max : m_dThisIterTopElevMax)                           \
    reduction(min : m_dThisIterTopElevMin)
#endif

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            m_ulThisIterNumCoastCells++;

         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         {
            // Is a sea cell
            m_dThisIterTotSeaDepth += m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();

            double dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
            m_pRasterGrid->m_Cell[nX][nY].IncrTotWaveAngle(dTmp);

            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
            m_pRasterGrid->m_Cell[nX][nY].IncrTotWaveHeight(dTmp);
         }

         double const dTopElev = m_pRasterGrid->m_Cell[nX][nY].dGetTopElevIncSea();

         // Get highest and lowest elevations of the top surface of the DEM
         if (dTopElev > m_dThisIterTopElevMax)
            m_dThisIterTopElevMax = dTopElev;

         if (dTopElev < m_dThisIterTopElevMin)
            m_dThisIterTopElevMin = dTopElev;

         // Reset the switch for platform erosion this timestep
         m_pRasterGrid->m_Cell[nX][nY].SetPlatformErosionThisIter(false);

         // First get fine talus depth for this cell, and allocate to a polygon total
         double dFineTalus = m_pRasterGrid->m_Cell[nX][nY].dGetFineTalusDepth();
         if (dFineTalus > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dFineTalus << " fine talus found on cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VdFineTalus[nCoast][nPoly] += dFineTalus;
         }

         // Next get sand talus depth for this cell, and allocate to a polygon total
         double dSandTalus = m_pRasterGrid->m_Cell[nX][nY].dGetSandTalusDepth();
         if (dSandTalus > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dSandTalus << " sand talus found on cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VdSandTalus[nCoast][nPoly] += dFineTalus;
         }

         // Finally get coarse talus depth for this cell, and allocate to a polygon total
         double dCoarseTalus = m_pRasterGrid->m_Cell[nX][nY].dGetCoarseTalusDepth();
         if (dCoarseTalus > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dCoarseTalus << " coarse talus found on cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VdSandTalus[nCoast][nPoly] += dFineTalus;
         }
      }
   }

   // No sea cells?
   if (m_ulThisIterNumSeaCells == 0)
      // All land, assume this is an error
      return RTN_ERR_NOSEACELLS;

   // Now go through all cells again and sort out suspended sediment load
   double const dSuspPerSeaCell = m_dThisIterFineSedimentToSuspension / static_cast<double>(m_ulThisIterNumSeaCells);

   // Parallelize the sediment distribution loop
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            m_pRasterGrid->m_Cell[nX][nY].AddSuspendedSediment(dSuspPerSeaCell);
      }
   }

   // Go along each coastline and update the grid with landform attributes, ready for next timestep
   for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
   {
      for (int j = 0; j < m_VCoast[i].nGetCoastlineSize(); j++)
      {
         int const nRet = nLandformToGrid(i, j);

         if (nRet != RTN_OK)
            return nRet;
      }
   }

   return RTN_OK;
}
