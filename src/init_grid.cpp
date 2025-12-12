/*!
   \file init_grid.cpp
   \brief Initialises the raster grid and calculates sea depth on each cell
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
#include <assert.h>

#include <climits>

#include <iostream>
using std::cerr;
using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cme.h"
#include "cell.h"
#include "coast.h"
#include "simulation.h"
#include "raster_grid.h"

//===============================================================================================================================
//! At the beginning of each timestep: clear vector coasts, profiles, and polygons, initialise the raster grid, and initialise some per-timestep accounting variables
//===============================================================================================================================
int CSimulation::nInitGridAndCalcStillWaterLevel(void)
{
   // Clear all vector coastlines, profiles, and polygons
   m_VCoast.clear();

   // m_VFloodWaveSetup.clear();
   m_VFloodWaveSetupSurge.clear();
   m_VFloodWaveSetupSurgeRunup.clear();

   // Do some every-timestep initialisation
   m_nXMinBoundingBox = INT_MAX;
   m_nXMaxBoundingBox = INT_MIN;
   m_nYMinBoundingBox = INT_MAX;
   m_nYMaxBoundingBox = INT_MIN;

   m_ulThisIterNumSeaCells = 0;
   m_ulThisIterNumCoastCells = 0;
   m_ulThisIterNumPotentialPlatformErosionCells = 0;
   m_ulThisIterNumActualPlatformErosionCells = 0;

   m_ulThisIterNumPotentialBeachErosionCells = 0;
   m_ulThisIterNumActualBeachErosionCells = 0;
   m_ulThisIterNumBeachDepositionCells = 0;

   m_dThisIterTotSeaDepth = 0;
   m_dThisIterPotentialPlatformErosion = 0;
   m_dThisIterPotentialBeachErosion = 0;
   m_dThisIterBeachErosionFine = 0;
   m_dThisIterBeachErosionSand = 0;
   m_dThisIterBeachErosionCoarse = 0;
   m_dThisIterBeachDepositionSand = 0;
   m_dThisIterBeachDepositionCoarse = 0;
   m_dThisIterPotentialSedLostBeachErosion = 0;
   m_dThisIterFineSedimentToSuspension = 0;
   m_dThisIterCliffCollapseErosionFineUncons = 0;
   m_dThisIterCliffCollapseErosionSandUncons = 0;
   m_dThisIterCliffCollapseErosionCoarseUncons = 0;
   m_dThisIterUnconsSandCliffDeposition = 0;
   m_dThisIterUnconsCoarseCliffDeposition = 0;
   m_dThisIterCliffCollapseErosionFineCons = 0;
   m_dThisIterCliffCollapseErosionSandCons = 0;
   m_dThisIterCliffCollapseErosionCoarseCons = 0;
   m_dThisIterActualPlatformErosionFineCons = 0;
   m_dThisIterActualPlatformErosionSandCons = 0;
   m_dThisIterActualPlatformErosionCoarseCons = 0;
   m_dThisIterLeftGridUnconsFine = 0; // TODO 067 Suspended fine sediment never decreases i.e. no suspended fine sediment ever leaves the grid. Is this OK?
   m_dThisIterLeftGridUnconsSand = 0;
   m_dThisIterLeftGridUnconsCoarse = 0;
   m_dThisiterUnconsFineInput = 0;
   m_dThisiterUnconsSandInput = 0;
   m_dThisiterUnconsCoarseInput = 0;

   for (int n = 0; n < m_nLayers; n++)
   {
      m_bConsChangedThisIter[n] = false;
      m_bUnconsChangedThisIter[n] = false;
   }

   // Re-calculate the depth of closure, in case deep water wave properties have changed
   CalcDepthOfClosure();

   int nZeroThickness = 0;

   m_dStartIterSuspFineAllCells = 0;
   m_dStartIterSuspFineInPolygons = 0;
   m_dStartIterUnconsFineAllCells = 0;
   m_dStartIterUnconsSandAllCells = 0;
   m_dStartIterUnconsCoarseAllCells = 0;
   m_dStartIterConsFineAllCells = 0;
   m_dStartIterConsSandAllCells = 0;
   m_dStartIterConsCoarseAllCells = 0;

   // And go through all cells in the RasterGrid array
   // Use OpenMP parallel loop with reduction clauses for thread-safe accumulation
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+ : nZeroThickness)                                            \
    reduction(+ : m_dStartIterConsFineAllCells, m_dStartIterConsSandAllCells, m_dStartIterConsCoarseAllCells) \
    reduction(+ : m_dStartIterSuspFineAllCells, m_dStartIterUnconsFineAllCells, m_dStartIterUnconsSandAllCells, m_dStartIterUnconsCoarseAllCells)
#endif

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         // Re-initialise values for this cell
         m_pRasterGrid->m_Cell[nX][nY].InitCell();

         if (m_ulIter == 1)
         {
            // For the first timestep only, check to see that all cells have some sediment on them
            double const dSedThickness = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedDepthAllLayers();

            if (dSedThickness <= 0)
            {
               nZeroThickness++;

               // Note: Logging from parallel regions can cause race conditions, but this is for debugging only. In production, consider collecting problematic cells and logging after the parallel region
               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               {
#ifdef _OPENMP
#pragma omp critical(logging)
#endif
                  LogStream << m_ulIter << ": " << WARN << "total sediment thickness is " << dSedThickness << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
               }
            }

            // For the first timestep only, calculate the elevation of all this cell's layers. During the rest of the simulation, each cell's elevation is re-calculated just after any change occurs on that cell
            m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
         }

         // Note that these totals include sediment which is both within and outside the polygons (because we have not yet defined polygons for this iteration, duh!)
         m_dStartIterConsFineAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetConsFineDepthAllLayers();
         m_dStartIterConsSandAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetConsSandDepthAllLayers();
         m_dStartIterConsCoarseAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetConsCoarseDepthAllLayers();

         m_dStartIterSuspFineAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();
         m_dStartIterUnconsFineAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetUnconsFineDepthAllLayers();
         m_dStartIterUnconsSandAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetUnconsSandDepthAllLayers();
         m_dStartIterUnconsCoarseAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetUnconsCoarseDepthAllLayers();

         if (m_bSingleDeepWaterWaveValues)
         {
            // If we have just a single measurement for deep water waves (either given by the user, or from a single wave station) then set all cells, even dry land cells, to the same value for deep water wave height, deep water wave orientation, and deep water period
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWaveHeight(m_dAllCellsDeepWaterWaveHeight);
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWaveAngle(m_dAllCellsDeepWaterWaveAngle);
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWavePeriod(m_dAllCellsDeepWaterWavePeriod);
         }
      }
   }

   if (m_bHaveWaveStationData && (! m_bSingleDeepWaterWaveValues))
   {
      // Each cell's value for deep water wave height and deep water wave orientation is interpolated from multiple user-supplied values
      int const nRet = nInterpolateAllDeepWaterWaveValues();

      if (nRet != RTN_OK)
         return nRet;

      /* for (int n = 0; n < m_VlDeepWaterWaveValuesAtTimestep.size(); n++)
         {
         if (m_ulIter == m_VlDeepWaterWaveValuesAtTimestep[n])
         {
            // OK, this timestep we are doing the calculation
            if (m_VlDeepWaterWaveValuesAtTimestep[n] > 1)
            {
               // TODO 036 For every timestep after the first, read in new values before doing the interpolation
            }

            // Interpolate values each cell's values for deep water height and orientation from user-supplied values
            int nRet = nInterpolateAllDeepWaterWaveValues();
            if (nRet != RTN_OK)
               return nRet;
         }
         }*/
   }

   if (nZeroThickness > 0)
   {
      cerr << m_ulIter << ": " << WARN << nZeroThickness << (nZeroThickness > 1 ? " cells" : " cell") << "( out of " << m_nXGridSize * m_nYGridSize << ") have no sediment, is this correct?" << endl;
      LogStream << m_ulIter << ": " << WARN << nZeroThickness << " cells have no sediment, is this correct?" << endl;
   }

   //    // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath;
   // strOutFile += "init_deep_water_wave_height_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_ulNumCells];
   // int nn = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   //          // Write this value to the array
   // pdRaster[nn] = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight();
   // nn++;
   // }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_nMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   //    // DEBUG CODE ===========================================================================================================

   //    // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "init_deep_water_wave_angle_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // nn = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   //          // Write this value to the array
   // pdRaster[nn] = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle();
   // nn++;
   // }
   // }
   //
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_nMissingValue);
   // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   //    // DEBUG CODE ===========================================================================================================

   //    // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "init_water_wave_angle_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // nn = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   //          // Write this value to the array
   // pdRaster[nn] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
   // nn++;
   // }
   // }
   //
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_nMissingValue);
   // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   //    // DEBUG CODE ===========================================================================================================

   //    // DEBUG CODE ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "init_water_wave_height_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // nn = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   //          // Write this value to the array
   // pdRaster[nn] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
   // nn++;
   // }
   // }
   //
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_nMissingValue);
   // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   //    // DEBUG CODE ===========================================================================================================

   return RTN_OK;
}
