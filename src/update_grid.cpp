/*!
   \file update_grid.cpp
   \brief Updates the raster grid
   \details TODO 001 A more detailed description of this routine.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <cfloat>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cme.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! Update all cells in the raster raster grid and do some per-timestep accounting
//===============================================================================================================================
int CSimulation::nUpdateGrid(void)
{
   // Go through all cells in the raster grid and calculate some this-timestep totals
   m_dThisIterTopElevMax = -DBL_MAX;
   m_dThisIterTopElevMin = DBL_MAX;

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
         }

         double const dTopElev = m_pRasterGrid->m_Cell[nX][nY].dGetTopElevIncSea();

         // Get highest and lowest elevations of the top surface of the DEM
         if (dTopElev > m_dThisIterTopElevMax)
            m_dThisIterTopElevMax = dTopElev;

         if (dTopElev < m_dThisIterTopElevMin)
            m_dThisIterTopElevMin = dTopElev;
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
