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

#include <vector>
using std::vector;

#include <cstddef>
using std::size_t;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cme.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

// #pragma omp declare reduction(matrix_add : vector<vector<double>> : \
//    for (size_t i = 0; i < omp_out.size(); ++i) \
//    { \
//       for (size_t j = 0; j < omp_out[i].size(); ++j) \
//       { \
//          omp_out[i][j] += omp_in[i][j]; \
//       } \
//    }) \
//    initializer(omp_priv = vector<vector<double>>(omp_orig.size(), vector<double>(omp_orig[0].size(), 0)))

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
   m_VVdFineTalus.resize(nCoastSize);
   m_VVdFineTalusAdded.resize(nCoastSize);
   m_VVdFineTalusRemoved.resize(nCoastSize);
   for (int n = 0; n < (nCoastSize); n++)
   {
      int nPoly = m_VCoast[n].nGetNumPolygons();
      m_VVdFineTalus[n].resize(nPoly);
      m_VVdFineTalusAdded[n].resize(nPoly);
      m_VVdFineTalusRemoved[n].resize(nPoly);
   }

   m_VVdSandTalus.resize(nCoastSize);
   m_VVdSandTalusAdded.resize(nCoastSize);
   m_VVdSandTalusRemoved.resize(nCoastSize);
   for (int n = 0; n < (nCoastSize); n++)
   {
      int nPoly = m_VCoast[n].nGetNumPolygons();
      m_VVdSandTalus[n].resize(nPoly);
      m_VVdSandTalusAdded[n].resize(nPoly);
      m_VVdSandTalusRemoved[n].resize(nPoly);
   }

   m_VVdCoarseTalus.resize(nCoastSize);
   m_VVdCoarseTalusAdded.resize(nCoastSize);
   m_VVdCoarseTalusRemoved.resize(nCoastSize);
   for (int n = 0; n < (nCoastSize); n++)
   {
      int nPoly = m_VCoast[n].nGetNumPolygons();
      m_VVdCoarseTalus[n].resize(nPoly);
      m_VVdCoarseTalusAdded[n].resize(nPoly);
      m_VVdCoarseTalusRemoved[n].resize(nPoly);
   }

   // Initialize reduction variables to zero
   m_ulThisIterNumCoastCells = 0;
   m_dThisIterTotSeaDepth = 0;

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2)                                 \
//     reduction(+ : m_ulThisIterNumCoastCells, m_dThisIterTotSeaDepth) \
//     reduction(max : m_dThisIterTopElevMax)                           \
//     reduction(min : m_dThisIterTopElevMin)
// #endif

   // #pragma omp parallel for reduction(matrix_add:m_VVdFineTalus)
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
               LogStream << m_ulIter << ":\t [" << nX << "][" << nY << "] fine talus = " << dFineTalus << " found but cell is not in a polygon" << endl;
            else
               m_VVdFineTalus[nCoast][nPoly] += dFineTalus;
         }

         // Get fine talus added to this cell, and allocate to a polygon total
         double dFineTalusAdded = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionFineToTalus();
         if (dFineTalusAdded > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dFineTalusAdded << " fine talus was added to cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VVdFineTalusAdded[nCoast][nPoly] += dFineTalusAdded;
         }

         // Get fine talus removed from this cell, and allocate to a polygon total
         double dFineTalusRemoved = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterFineTalusToSuspension();
         if (dFineTalusRemoved > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dFineTalusRemoved << " fine talus was removed from cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VVdFineTalusRemoved[nCoast][nPoly] += dFineTalusRemoved;
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
               m_VVdSandTalus[nCoast][nPoly] += dSandTalus;
         }

         // Get sand talus added to this cell, and allocate to a polygon total
         double dSandTalusAdded = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionSandToTalus();
         if (dSandTalusAdded > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dSandTalusAdded << " sand talus was added to cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VVdSandTalusAdded[nCoast][nPoly] += dSandTalusAdded;
         }

         // Get sand talus removed from this cell, and allocate to a polygon total
         double dSandTalusRemoved = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterSandTalusToUncons();
         if (dSandTalusRemoved > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dSandTalusRemoved << " sand talus was removed from cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VVdSandTalusRemoved[nCoast][nPoly] += dSandTalusRemoved;
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
               m_VVdCoarseTalus[nCoast][nPoly] += dCoarseTalus;
         }

         // Get coarse talus added to this cell, and allocate to a polygon total
         double dCoarseTalusAdded = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionCoarseToTalus();
         if (dCoarseTalusAdded > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dCoarseTalusAdded << " coarse talus was added to cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VVdCoarseTalusAdded[nCoast][nPoly] += dCoarseTalusAdded;
         }

         // Get coarse talus removed from this cell, and allocate to a polygon total
         double dCoarseTalusRemoved = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCoarseTalusToUncons();
         if (dCoarseTalusRemoved > 0)
         {
            int nCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
            int nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            if ((nPoly == INT_NODATA) || (nCoast == INT_NODATA))
               LogStream << m_ulIter << ":\t " << dCoarseTalusRemoved << " coarse talus was removed from cell [" << nX << "][" << nY << "] but cell is not in a polygon" << endl;
            else
               m_VVdCoarseTalusRemoved[nCoast][nPoly] += dCoarseTalusRemoved;
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
// #ifdef _OPENMP
// #pragma omp parallel for collapse(2)
// #endif

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
