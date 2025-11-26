/*!
   \file gis_raster.cpp
   \brief These functions use GDAL (at least version 2) to read and write raster
   GIS files in several formats
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
*/

/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <assert.h>

#include <cstdio>

#include <cmath>
using std::atan2;
using std::hypot;
using std::isfinite;
using std::isnan;
using std::sqrt;

#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::stringstream;

#include <string>
using std::to_string;

#include <cpl_conv.h>
#include <cpl_error.h>
#include <cpl_string.h>
#include <gdal.h>
#include <gdal_alg.h>
#include <gdal_priv.h>

#include "2di_point.h"
#include "cme.h"
#include "coast.h"
#include "simulation.h"
#include "spatial_interpolation.h"

//===============================================================================================================================
//! Initialize GDAL with performance optimizations
//===============================================================================================================================
void CSimulation::InitializeGDALPerformance(void)
{
   // Enable GDAL threading - use all available CPU cores
#ifdef _OPENMP
   CPLSetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS");
#else
   CPLSetConfigOption("GDAL_NUM_THREADS", "4");                            // Fallback for non-OpenMP builds
#endif

   // Optimize GDAL memory usage and caching
   CPLSetConfigOption("GDAL_CACHEMAX", "1.5GB");                           // 2GB cache for large grids (was 1GB)
   CPLSetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR");        // Faster file access
   CPLSetConfigOption("VSI_CACHE", "TRUE");                                // Enable virtual file system cache
   CPLSetConfigOption("VSI_CACHE_SIZE", "512MB");                          // Was 256MB VSI cache

   // Block and chunk optimizations for raster operations
   CPLSetConfigOption("GDAL_TIFF_INTERNAL_MASK_TO_8BIT", "YES");
   CPLSetConfigOption("GDAL_RASTERIO_RESAMPLING", "CUBIC");                // Better for coastal DEM data

   // Grid creation optimizations (for GDALGridCreate performance)
   CPLSetConfigOption("GDAL_GRID_MAX_POINTS_PER_QUADTREE_LEAF", "1024");   // Increased from 512
   CPLSetConfigOption("GDAL_GRID_POINT_COUNT_THRESHOLD", "100");           // New 2024 option

   // Thread-safe dataset access (GDAL 3.10+)
   CPLSetConfigOption("GDAL_DATASET_CACHE_SIZE", "64");                    // Cache more datasets

   // Compression optimizations for output
   CPLSetConfigOption("GDAL_TIFF_OVR_BLOCKSIZE", "512");                   // Optimal for coastal data

   // Memory allocator optimization for multi-threading
   CPLSetConfigOption("CPL_VSIL_USE_TEMP_FILE_FOR_RANDOM_WRITE", "YES");

   // Disable GDAL warnings for cleaner output (optional)
   // CPLSetConfigOption("CPL_LOG", "/dev/null");
   // Debugging (remove in production)
   // CPLSetConfigOption("CPL_DEBUG", "ON");
   // CPLSetConfigOption("GDAL_DEBUG", "ON");

   m_bGDALOptimisations = true;
}

//===============================================================================================================================
//! Reads a raster DEM of basement elevation data to the Cell array
//===============================================================================================================================
int CSimulation::nReadRasterBasementDEM(void)
{
   // Initialize GDAL performance settings (only needs to be done once)
   static bool bGDALInitialized = false;

   if (! bGDALInitialized)
   {
      InitializeGDALPerformance();
      bGDALInitialized = true;
   }

   // Use GDAL to create a dataset object, which then opens the DEM file
   GDALDataset *pGDALDataset = static_cast<GDALDataset *>(GDALOpen(m_strInitialBasementDEMFile.c_str(), GA_ReadOnly));

   if (NULL == pGDALDataset)
   {
      // Can't open file (note will already have sent GDAL error message to stdout)
      cerr << ERR << "cannot open " << m_strInitialBasementDEMFile << " for input: " << CPLGetLastErrorMsg() << endl;
      return RTN_ERR_DEMFILE;
   }

   // Opened OK, so get GDAL basement DEM dataset information
   m_strGDALBasementDEMDriverCode = pGDALDataset->GetDriver()->GetDescription();
   m_strGDALBasementDEMDriverDesc = pGDALDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME);
   m_strGDALBasementDEMProjection = pGDALDataset->GetProjectionRef();

   if (m_strGDALBasementDEMProjection.empty())
   {
      // TODO
      m_strGDALBasementDEMProjection = "";
      // pGDALDataset->SetProjectionRef(m_strGDALBasementDEMProjection);

      // ENGCRS["Plane",
      //     EDATUM["Unknown engineering datum"],
      //     CS[Cartesian,2],
      //         AXIS["(E)",east,
      //             ORDER[1],
      //             LENGTHUNIT["Meter",1]],
      //         AXIS["(N)",north,
      //             ORDER[2],
      //             LENGTHUNIT["Meter",1]]]
   }
   else
   {
      // We have reference units, so check that they are in metres (note US spelling)
      if (! m_strGDALBasementDEMProjection.empty())
      {
         string const strTmp = strToLower(&m_strGDALBasementDEMProjection);

         if ((strTmp.find("meter") == string::npos) && (strTmp.find("metre") == string::npos))
         {
            // error: x-y values must be in metres
            cerr << ERR << "GIS file x-y values (" << m_strGDALBasementDEMProjection << ") in " << m_strInitialBasementDEMFile << " must be in metres" << endl;
            return RTN_ERR_DEMFILE;
         }
      }
   }

   // Now get dataset size, and do some rudimentary checks
   m_nXGridSize = pGDALDataset->GetRasterXSize();

   if (m_nXGridSize == 0)
   {
      // Error: silly number of columns specified
      cerr << ERR << "invalid number of columns (" << m_nXGridSize << ") in " << m_strInitialBasementDEMFile << endl;
      return RTN_ERR_DEMFILE;
   }

   m_nYGridSize = pGDALDataset->GetRasterYSize();

   if (m_nYGridSize == 0)
   {
      // Error: silly number of rows specified
      cerr << ERR << "invalid number of rows (" << m_nYGridSize << ") in " << m_strInitialBasementDEMFile << endl;
      return RTN_ERR_DEMFILE;
   }

   // Get geotransformation info (see http://www.gdal.org/classGDALDataset.html)
   if (CE_Failure == pGDALDataset->GetGeoTransform(m_dGeoTransform))
   {
      // Can't get geotransformation (note will already have sent GDAL error message to stdout)
      cerr << ERR << CPLGetLastErrorMsg() << " in " << m_strInitialBasementDEMFile << endl;
      return RTN_ERR_DEMFILE;
   }

   // CoastalME can only handle rasters that are oriented N-S and W-E. (If you need to work with a raster that is oriented differently, then you must rotate it before running CoastalME). So here we check whether row rotation (m_dGeoTransform[2]) and column rotation (m_dGeoTransform[4]) are both zero. See https://gdal.org/tutorials/geotransforms_tut.html
   if ((! bFPIsEqual(m_dGeoTransform[2], 0.0, TOLERANCE)) || (! bFPIsEqual(m_dGeoTransform[4], 0.0, TOLERANCE)))
   {
      // Error: not oriented NS and W-E
      cerr << ERR << m_strInitialBasementDEMFile << " is not oriented N-S and W-E. Row rotation = " << m_dGeoTransform[2] << " and column rotation = " << m_dGeoTransform[4] << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   // Get the X and Y cell sizes, in external CRS units. Note that while the cell is supposed to be square, it may not be exactly so due to oddities with some GIS calculations
   double const dCellSideX = tAbs(m_dGeoTransform[1]);
   double const dCellSideY = tAbs(m_dGeoTransform[5]);

   // Check that the cell is more or less square
   if (! bFPIsEqual(dCellSideX, dCellSideY, 1e-2))
   {
      // Error: cell is not square enough
      cerr << ERR << "cell is not square in " << m_strInitialBasementDEMFile << ", is " << dCellSideX << " x " << dCellSideY << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   // Calculate the average length of cell side, the cell's diagonal, and the area of a cell (in external CRS units)
   m_dCellSide = (dCellSideX + dCellSideY) / 2.0;
   m_dCellArea = m_dCellSide * m_dCellSide;
   m_dCellDiagonal = hypot(m_dCellSide, m_dCellSide);

   // And calculate the inverse values
   m_dInvCellSide = 1 / m_dCellSide;
   m_dInvCellDiagonal = 1 / m_dCellDiagonal;

   // Save some values in external CRS
   m_dNorthWestXExtCRS = m_dGeoTransform[0] - (m_dGeoTransform[1] / 2);
   m_dNorthWestYExtCRS = m_dGeoTransform[3] - (m_dGeoTransform[5] / 2);
   m_dSouthEastXExtCRS = m_dGeoTransform[0] + (m_nXGridSize * m_dGeoTransform[1]) + (m_dGeoTransform[1] / 2);
   m_dSouthEastYExtCRS = m_dGeoTransform[3] + (m_nYGridSize * m_dGeoTransform[5]) + (m_dGeoTransform[5] / 2);

   // And calc the grid area in external CRS units
   m_dExtCRSGridArea = tAbs(m_dNorthWestXExtCRS - m_dSouthEastXExtCRS) * tAbs(m_dNorthWestYExtCRS * m_dSouthEastYExtCRS);

   // Now get GDAL raster band information
   GDALRasterBand *pGDALBand = pGDALDataset->GetRasterBand(1);
   int nBlockXSize = 0, nBlockYSize = 0;
   pGDALBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
   m_strGDALBasementDEMDataType = GDALGetDataTypeName(pGDALBand->GetRasterDataType());

   // If we have value units, then check them
   string const strUnits = pGDALBand->GetUnitType();

   if ((! strUnits.empty()) && (strUnits.find('m') == string::npos))
   {
      // Error: value units must be m
      cerr << ERR << "DEM vertical units are (" << strUnits << " ) in " << m_strInitialBasementDEMFile << ", should be 'm'" << endl;
      return RTN_ERR_DEMFILE;
   }

   // If present, get the missing value (NODATA) setting
   CPLPushErrorHandler(CPLQuietErrorHandler);                     // Needed to get next line to fail silently, if it fails
   double const dMissingValue = pGDALBand->GetNoDataValue();      // Will fail for some formats
   CPLPopErrorHandler();

   if (! bFPIsEqual(dMissingValue, m_dMissingValue, TOLERANCE))
   {
      cerr << NOTE << "NODATA value in " << m_strInitialBasementDEMFile << " is " << dMissingValue << "\n         instead using CoastalME's default floating-point NODATA value " << m_dMissingValue << endl;
   }

   // Next allocate memory for a 2D array of raster cell objects: tell the user what is happening
   AnnounceAllocateMemory();
   int const nRet = m_pRasterGrid->nCreateGrid();

   if (nRet != RTN_OK)
      return nRet;

   // Allocate memory for a 1D floating-point array, to hold the scan line for GDAL
   double *pdScanline = new double[m_nXGridSize];

   if (NULL == pdScanline)
   {
      // Error, can't allocate memory
      cerr << ERR << "cannot allocate memory for " << m_nXGridSize << " x 1D array" << endl;
      return (RTN_ERR_MEMALLOC);
   }

   // Now read in the data
   for (int j = 0; j < m_nYGridSize; j++)
   {
      // Read scanline
      if (CE_Failure == pGDALBand->RasterIO(GF_Read, 0, j, m_nXGridSize, 1, pdScanline, m_nXGridSize, 1, GDT_Float64, 0, 0, NULL))
      {
         // Error while reading scanline
         cerr << ERR << CPLGetLastErrorMsg() << " in " << m_strInitialBasementDEMFile << endl;
         return RTN_ERR_DEMFILE;
      }

      // All OK, so read scanline into cell elevations (including any missing values)
      for (int i = 0; i < m_nXGridSize; i++)
      {
         double dTmp = pdScanline[i];

         if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
            dTmp = m_dMissingValue;

         m_pRasterGrid->m_Cell[i][j].SetBasementElev(dTmp);
      }
   }

   // Finished, so get rid of dataset object
   GDALClose(pGDALDataset);

   // Get rid of memory allocated to this array
   delete[] pdScanline;

   return RTN_OK;
}

//===============================================================================================================================
//! Mark cells which are at the edge of a bounding box which represents the valid part of the grid, as defined by the basement layer. The valid part of the grid may be the whole grid, or only part of the whole grid. The bounding box may be an irregular shape (but may not have re-entrant edges): simple shapes are more likely to work correctly
//===============================================================================================================================
int CSimulation::nMarkBoundingBoxEdgeCells(void)
{
   // The bounding box must touch the edge of the grid at least once on each side of the grid. Store these corner points in the sequence [NW NE SW SE]. Start with the top (north) edge, search NW to NE (left to right)
   bool bFound = false;
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      if (bFound)
         break;

      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            // This is the NW corner
            CGeom2DIPoint const PtiTmp(nX, nY);
            m_VPtiBoundingBoxCorner.push_back(PtiTmp);

            bFound = true;
            break;
         }
      }
   }

   if (! bFound)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": NW (top left) corner of bounding box not found" << endl;

      return RTN_ERR_BOUNDING_BOX;
   }

   // Again go along the top (north) edge, but this time search NE to NW (right to left)
   bFound = false;
   for (int nX = m_nXGridSize-1; nX >= 0; nX--)
   {
      if (bFound)
         break;

      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            // This is the NE corner
            CGeom2DIPoint const PtiTmp(nX, nY);
            m_VPtiBoundingBoxCorner.push_back(PtiTmp);

            bFound = true;
            break;
         }
      }
   }

   if (! bFound)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": NE (top right) edge of bounding box not found" << endl;

      return RTN_ERR_BOUNDING_BOX;
   }


   // Now go along the bottom (south) edge, search SW to SE (top to bottom)
   bFound = false;
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      if (bFound)
         break;

      for (int nY = m_nYGridSize-1; nY >= 0; nY--)
      {
         if (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            // This is the SW corner
            CGeom2DIPoint const PtiTmp(nX, nY);
            m_VPtiBoundingBoxCorner.push_back(PtiTmp);

            bFound = true;
            break;
         }
      }
   }

   if (! bFound)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": SW (bottom left) corner of bounding box not found" << endl;

      return RTN_ERR_BOUNDING_BOX;
   }

   // Finally, again go along the bottom (south) edge, but this time search SE to SW (bottom to top)
   bFound = false;
   for (int nX = m_nXGridSize-1; nX >= 0; nX--)
   {
      if (bFound)
         break;

      for (int nY = m_nYGridSize-1; nY >= 0; nY--)
      {
         if (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            // This is the NE corner
            CGeom2DIPoint const PtiTmp(nX, nY);
            m_VPtiBoundingBoxCorner.push_back(PtiTmp);

            bFound = true;
            break;
         }
      }
   }

   if (! bFound)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": SE (bottom right) edge of bounding box not found" << endl;

      return RTN_ERR_BOUNDING_BOX;
   }

   // OK, so we have a point at each corner [NW NE SW SE] of the grid, so start at this point and find the edges of the bounding box
   for (int nX = m_VPtiBoundingBoxCorner[0].nGetX(); nX <= m_VPtiBoundingBoxCorner[1].nGetX(); nX++)
   {
      // Search left to right (W to E)
      bFound = false;
      for (int nY = m_VPtiBoundingBoxCorner[0].nGetY(); nY < m_nYGridSize; nY++)
      {
         // Search top to bottom (N to S)
         if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            m_ulMissingValueBasementCells++;
            continue;
         }

         // Found a bounding box edge cell
         m_pRasterGrid->m_Cell[nX][nY].SetBoundingBoxEdge(NORTH);
         m_VPtiNorthEdgeCell.push_back(CGeom2DIPoint(nX, nY));

         bFound = true;
         break;
      }

      if (! bFound)
      {
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": could not find a bounding box edge cell for grid column " << nX << endl;

         return RTN_ERR_BOUNDING_BOX;
      }
   }

   // Right (east) edge
   for (int nY = m_VPtiBoundingBoxCorner[1].nGetY(); nY <= m_VPtiBoundingBoxCorner[3].nGetY(); nY++)
   {
      // Search top to bottom (N to S)
      bFound = false;
      for (int nX = m_VPtiBoundingBoxCorner[1].nGetX(); nX >= 0; nX--)
      {
         // Search right to left (E to W)
         if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            m_ulMissingValueBasementCells++;
            continue;
         }

         // Found a bounding box edge cell
         m_pRasterGrid->m_Cell[nX][nY].SetBoundingBoxEdge(EAST);
         m_VPtiEastEdgeCell.push_back(CGeom2DIPoint(nX, nY));

         bFound = true;
         break;
      }

      if (! bFound)
      {
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": could not find a bounding box edge cell for grid row " << nY << endl;

         return RTN_ERR_BOUNDING_BOX;
      }
   }

   // Bottom (south) edge
   for (int nX = m_VPtiBoundingBoxCorner[2].nGetX(); nX >= m_VPtiBoundingBoxCorner[3].nGetX(); nX++)
   {
      // Searching left to right (W to E)
      bFound = false;
      for (int nY = m_VPtiBoundingBoxCorner[2].nGetY(); nY >= 0; nY--)
      {
         // Searching bottom to top (S to N)
         if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            m_ulMissingValueBasementCells++;
            continue;
         }

         // Found a bounding box edge cell
         m_pRasterGrid->m_Cell[nX][nY].SetBoundingBoxEdge(SOUTH);
         m_VPtiSouthEdgeCell.push_back(CGeom2DIPoint(nX, nY));

         bFound = true;
         break;
      }

      if (! bFound)
      {
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": could not find a bounding box edge cell for grid column " << nX << endl;

         return RTN_ERR_BOUNDING_BOX;
      }
   }

   // Left (west) edge
   for (int nY = m_VPtiBoundingBoxCorner[0].nGetY(); nY >= m_VPtiBoundingBoxCorner[2].nGetY(); nY++)
   {
      // Search top to bottom (N to S)
      bFound = false;
      for (int nX = m_VPtiBoundingBoxCorner[0].nGetX(); nX < m_nXGridSize - 1; nX++)
      {
         // Searching left to right (W to E)
         if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            m_ulMissingValueBasementCells++;
            continue;
         }

         // Found a bounding box edge cell
         m_pRasterGrid->m_Cell[nX][nY].SetBoundingBoxEdge(WEST);
         m_VPtiWestEdgeCell.push_back(CGeom2DIPoint(nX, nY));

         bFound = true;
         break;
      }

      if (! bFound)
      {
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": could not find a bounding box edge cell for grid row " << nY << endl;

         return RTN_ERR_BOUNDING_BOX;
      }
   }

   // Finally concatenate the four edge cell vectors to make an all-edge-cells vector, in clockwise sequence
   m_VPtiAllEdgeCell.insert(m_VPtiAllEdgeCell.end(), m_VPtiNorthEdgeCell.begin(), m_VPtiNorthEdgeCell.end());
   m_VPtiAllEdgeCell.insert(m_VPtiAllEdgeCell.end(), m_VPtiEastEdgeCell.begin(), m_VPtiEastEdgeCell.end());
   m_VPtiAllEdgeCell.insert(m_VPtiAllEdgeCell.end(), m_VPtiSouthEdgeCell.end(), m_VPtiSouthEdgeCell.begin());     // Reverse
   m_VPtiAllEdgeCell.insert(m_VPtiAllEdgeCell.end(), m_VPtiWestEdgeCell.end(), m_VPtiWestEdgeCell.begin());       // Reverse

   return RTN_OK;
}

//===============================================================================================================================
//! Find the elevations of the four bounding-box corners, and use these to choose a sea-to-land (i.e. low elevation to high elevation) direction of processing for each bounding-box edge
//===============================================================================================================================
void CSimulation::CalcGridEdgeSeaToLandDirection(void)
{
   // The bounding box corner points are in the sequence [NW NE SW SE]
   vector<double> VdBoundingBoxCornerElev;

   for (int i = 0; i < 4; i++)
   {
      int nX = m_VPtiBoundingBoxCorner[i].nGetX();
      int nY = m_VPtiBoundingBoxCorner[i].nGetY();
      m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
      VdBoundingBoxCornerElev.push_back(m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus());
   }

   // We have the elevations of the four bounding box corner points, in the sequence [NW NE SW SE]. When finding possible start/finish points for coast tracing later on, we wish to process each side of the grid in the sequence low elevation (sea) to high elevation (land). So set these directions now
   if (VdBoundingBoxCornerElev[0] > VdBoundingBoxCornerElev[1])
   {
      // NW (left) point is higher than NE (right) point, so search north edge for coast points right to left
      m_bSearchNorthEdgeForward = false;
   }
   else
   {
      // NW (left) point is lower than NE (right) point, so search north edge for coast points left to right
      m_bSearchNorthEdgeForward = true;
   }

   if (VdBoundingBoxCornerElev[1] > VdBoundingBoxCornerElev[3])
   {
      // NE (top) point is higher than SE (bottom) point, so search east edge for coast points bottom to top
      m_bSearchEastEdgeForward = false;
   }
   else
   {
      // NE (bottom) point is lower than SE (top) point, so search east edge for coast points top to bottom
      m_bSearchEastEdgeForward = true;
   }

   if (VdBoundingBoxCornerElev[2] > VdBoundingBoxCornerElev[3])
   {
      // SW (left) point is higher than SE (right) point, so search south edge for coast points right to left
      m_bSearchSouthEdgeForward = false;
   }
   else
   {
      // SW (left) point is lower than SE (right) point, so search south edge for coast points left to right
      m_bSearchSouthEdgeForward = true;
   }

   if (VdBoundingBoxCornerElev[0] > VdBoundingBoxCornerElev[2])
   {
      // NW (top) point is higher than SW (bottom) point, so search west edge for coast points bottom to top
      m_bSearchWestEdgeForward = false;
   }
   else
   {
      // NW (top) point is lower than SW (bottom) point, so search west edge for coast points top to bottom
      m_bSearchWestEdgeForward = true;
   }
}

//===============================================================================================================================
//! Reads raster GIS datafiles into the RasterGrid array
//===============================================================================================================================
int CSimulation::nReadRasterGISFile(int const nDataItem, int const nLayer)
{
   string strGISFile;
   string strDriverCode;
   string strDriverDesc;
   string strProjection;
   string strDataType;

   switch (nDataItem)
   {
      case (LANDFORM_RASTER):
         // Initial Landform Class GIS data
         strGISFile = m_strInitialLandformFile;
         break;

      case (INTERVENTION_CLASS_RASTER):
         // Intervention class
         strGISFile = m_strInterventionClassFile;
         break;

      case (INTERVENTION_HEIGHT_RASTER):
         // Intervention height
         strGISFile = m_strInterventionHeightFile;
         break;

      case (SUSP_SED_RASTER):
         // Initial Suspended Sediment GIS data
         strGISFile = m_strInitialSuspSedimentFile;
         break;

      case (FINE_UNCONS_RASTER):
         // Initial Unconsolidated Fine Sediment GIS data
         strGISFile = m_VstrInitialFineUnconsSedimentFile[nLayer];
         break;

      case (SAND_UNCONS_RASTER):
         // Initial Unconsolidated Sand Sediment GIS data
         strGISFile = m_VstrInitialSandUnconsSedimentFile[nLayer];
         break;

      case (COARSE_UNCONS_RASTER):
         // Initial Unconsolidated Coarse Sediment GIS data
         strGISFile = m_VstrInitialCoarseUnconsSedimentFile[nLayer];
         break;

      case (FINE_CONS_RASTER):
         // Initial Consolidated Fine Sediment GIS data
         strGISFile = m_VstrInitialFineConsSedimentFile[nLayer];
         break;

      case (SAND_CONS_RASTER):
         // Initial Consolidated Sand Sediment GIS data
         strGISFile = m_VstrInitialSandConsSedimentFile[nLayer];
         break;

      case (COARSE_CONS_RASTER):
         // Initial Consolidated Coarse Sediment GIS data
         strGISFile = m_VstrInitialCoarseConsSedimentFile[nLayer];
         break;
   }

   // Use GDAL to create a dataset object, which then opens the GIS file
   GDALDataset *pGDALDataset =
      static_cast<GDALDataset *>(GDALOpen(strGISFile.c_str(), GA_ReadOnly));

   if (NULL == pGDALDataset)
   {
      // Can't open file (note will already have sent GDAL error message to
      // stdout)
      cerr << ERR << "cannot open " << strGISFile
           << " for input: " << CPLGetLastErrorMsg() << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   // Opened OK, so get dataset information
   strDriverCode = pGDALDataset->GetDriver()->GetDescription();
   strDriverDesc = pGDALDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME);
   strProjection = pGDALDataset->GetProjectionRef();

   // Get geotransformation info
   double dGeoTransform[6];
   if (CE_Failure == pGDALDataset->GetGeoTransform(dGeoTransform))
   {
      // Can't get geotransformation (note will already have sent GDAL error
      // message to stdout)
      cerr << ERR << CPLGetLastErrorMsg() << " in " << strGISFile << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   // Now get dataset size, and do some checks
   int const nTmpXSize = pGDALDataset->GetRasterXSize();
   if (nTmpXSize != m_nXGridSize)
   {
      // Error: incorrect number of columns specified
      cerr << ERR << "different number of columns in " << strGISFile << " (" << nTmpXSize << ") and " << m_strInitialBasementDEMFile << "(" << m_nXGridSize << ")" << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   int const nTmpYSize = pGDALDataset->GetRasterYSize();
   if (nTmpYSize != m_nYGridSize)
   {
      // Error: incorrect number of rows specified
      cerr << ERR << "different number of rows in " << strGISFile << " (" << nTmpYSize << ") and " << m_strInitialBasementDEMFile << " (" << m_nYGridSize << ")" << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   double dTmp = m_dGeoTransform[0] - (m_dGeoTransform[1] / 2);
   if (! bFPIsEqual(dTmp, m_dNorthWestXExtCRS, TOLERANCE))
   {
      // Error: different min x from DEM file
      cerr << ERR << "different min x values in " << strGISFile << " (" << dTmp << ") and " << m_strInitialBasementDEMFile << " (" << m_dNorthWestXExtCRS << ")" << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   dTmp = m_dGeoTransform[3] - (m_dGeoTransform[5] / 2);
   if (! bFPIsEqual(dTmp, m_dNorthWestYExtCRS, TOLERANCE))
   {
      // Error: different min x from DEM file
      cerr << ERR << "different min y values in " << strGISFile << " (" << dTmp << ") and " << m_strInitialBasementDEMFile << " ("  << m_dNorthWestYExtCRS << ")" << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   double const dTmpResX = tAbs(dGeoTransform[1]);
   if (! bFPIsEqual(dTmpResX, m_dCellSide, 1e-2))
   {
      // Error: different cell size in X direction: note that due to rounding errors in some GIS packages, must expect some discrepancies
      cerr << ERR << "cell size in X direction (" << dTmpResX << ") in " << strGISFile << " differs from cell size in of basement DEM (" << m_dCellSide << ")" << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   double const dTmpResY = tAbs(dGeoTransform[5]);
   if (! bFPIsEqual(dTmpResY, m_dCellSide, 1e-2))
   {
      // Error: different cell size in Y direction: note that due to rounding errors in some GIS packages, must expect some discrepancies
      cerr << ERR << "cell size in Y direction (" << dTmpResY << ") in " << strGISFile << " differs from cell size of basement DEM (" << m_dCellSide << ")" << endl;
      return (RTN_ERR_RASTER_FILE_READ);
   }

   // Now get GDAL raster band information
   GDALRasterBand *pGDALBand = pGDALDataset->GetRasterBand(1);      // TODO 028 Give a message if there are several bands
   int nBlockXSize = 0, nBlockYSize = 0;
   pGDALBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
   strDataType = GDALGetDataTypeName(pGDALBand->GetRasterDataType());

   switch (nDataItem)
   {
      case (LANDFORM_RASTER):
         // Initial Landform Class GIS data
         m_strGDALLDriverCode = strDriverCode;
         m_strGDALLDriverDesc = strDriverDesc;
         m_strGDALLProjection = strProjection;
         m_strGDALLDataType = strDataType;
         break;

      case (INTERVENTION_CLASS_RASTER):
         // Intervention class
         m_strGDALICDriverCode = strDriverCode;
         m_strGDALICDriverDesc = strDriverDesc;
         m_strGDALICProjection = strProjection;
         m_strGDALICDataType = strDataType;
         break;

      case (INTERVENTION_HEIGHT_RASTER):
         // Intervention height
         m_strGDALIHDriverCode = strDriverCode;
         m_strGDALIHDriverDesc = strDriverDesc;
         m_strGDALIHProjection = strProjection;
         m_strGDALIHDataType = strDataType;
         break;

      case (SUSP_SED_RASTER):
         // Initial Suspended Sediment GIS data
         m_strGDALISSDriverCode = strDriverCode;
         m_strGDALISSDriverDesc = strDriverDesc;
         m_strGDALISSProjection = strProjection;
         m_strGDALISSDataType = strDataType;
         break;

      case (FINE_UNCONS_RASTER):
         // Initial Unconsolidated Fine Sediment GIS data
         m_VstrGDALIUFDriverCode[nLayer] = strDriverCode;
         m_VstrGDALIUFDriverDesc[nLayer] = strDriverDesc;
         m_VstrGDALIUFProjection[nLayer] = strProjection;
         m_VstrGDALIUFDataType[nLayer] = strDataType;
         break;

      case (SAND_UNCONS_RASTER):
         // Initial Unconsolidated Sand Sediment GIS data
         m_VstrGDALIUSDriverCode[nLayer] = strDriverCode;
         m_VstrGDALIUSDriverDesc[nLayer] = strDriverDesc;
         m_VstrGDALIUSProjection[nLayer] = strProjection;
         m_VstrGDALIUSDataType[nLayer] = strDataType;
         break;

      case (COARSE_UNCONS_RASTER):
         // Initial Unconsolidated Coarse Sediment GIS data
         m_VstrGDALIUCDriverCode[nLayer] = strDriverCode;
         m_VstrGDALIUCDriverDesc[nLayer] = strDriverDesc;
         m_VstrGDALIUCProjection[nLayer] = strProjection;
         m_VstrGDALIUCDataType[nLayer] = strDataType;
         break;

      case (FINE_CONS_RASTER):
         // Initial Consolidated Fine Sediment GIS data
         m_VstrGDALICFDriverCode[nLayer] = strDriverCode;
         m_VstrGDALICFDriverDesc[nLayer] = strDriverDesc;
         m_VstrGDALICFProjection[nLayer] = strProjection;
         m_VstrGDALICFDataType[nLayer] = strDataType;
         break;

      case (SAND_CONS_RASTER):
         // Initial Consolidated Sand Sediment GIS data
         m_VstrGDALICSDriverCode[nLayer] = strDriverCode;
         m_VstrGDALICSDriverDesc[nLayer] = strDriverDesc;
         m_VstrGDALICSProjection[nLayer] = strProjection;
         m_VstrGDALICSDataType[nLayer] = strDataType;
         break;

      case (COARSE_CONS_RASTER):
         // Initial Consolidated Coarse Sediment GIS data
         m_VstrGDALICCDriverCode[nLayer] = strDriverCode;
         m_VstrGDALICCDriverDesc[nLayer] = strDriverDesc;
         m_VstrGDALICCProjection[nLayer] = strProjection;
         m_VstrGDALICCDataType[nLayer] = strDataType;
         break;
   }

   // If present, get the missing value setting
   string const strTmp = strToLower(&strDataType);
   if (strTmp.find("int") != string::npos)
   {
      // This is an integer layer
      CPLPushErrorHandler(CPLQuietErrorHandler);                                 // Needed to get next line to fail silently, if it fails
      m_nGISMissingValue = static_cast<int>(pGDALBand->GetNoDataValue());        // Note will fail for some formats
      CPLPopErrorHandler();

      if (m_nGISMissingValue != m_nMissingValue)
      {
         cerr << NOTE << "NODATA value in " << strGISFile << " is " << m_nGISMissingValue << "\n         instead using CoatalME's default integer NODATA value " << m_nMissingValue << endl;
      }
   }
   else
   {
      // This is a floating point layer
      CPLPushErrorHandler(CPLQuietErrorHandler);                                 // Needed to get next line to fail silently, if it fails
      m_dGISMissingValue = pGDALBand->GetNoDataValue();                          // Note will fail for some formats
      CPLPopErrorHandler();

      if (! bFPIsEqual(m_dGISMissingValue, m_dMissingValue, TOLERANCE))
      {
         cerr << NOTE << "NODATA value in " << strGISFile << " is " << m_dGISMissingValue << "\n         instead using CoastalME's default floating-point NODATA value " << m_dMissingValue << endl;
      }
   }

   // Allocate memory for a 1D array, to hold the scan line for GDAL
   double *pdScanline = new double[m_nXGridSize];
   if (NULL == pdScanline)
   {
      // Error, can't allocate memory
      cerr << ERR << "cannot allocate memory for " << m_nXGridSize << " x 1D array" << endl;
      return (RTN_ERR_MEMALLOC);
   }

   // Now read in the data
   int nMissing = 0;

   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      // Read scanline
      if (CE_Failure == pGDALBand->RasterIO(GF_Read, 0, nY, m_nXGridSize, 1, pdScanline, m_nXGridSize, 1, GDT_Float64, 0, 0, NULL))
      {
         // Error while reading scanline
         cerr << ERR << CPLGetLastErrorMsg() << " in " << strGISFile << endl;
         return (RTN_ERR_RASTER_FILE_READ);
      }

      // All OK, so read scanline into cells (including any missing values)
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         int nTmp;

         switch (nDataItem)
         {
            case (LANDFORM_RASTER):
               // Initial Landform Class GIS data, is integer TODO 030 Do we also need
               // a landform sub-category input?
               nTmp = static_cast<int>(pdScanline[nX]);

               if ((isnan(nTmp)) || (nTmp == m_nGISMissingValue))
               {
                  nTmp = m_nMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(nTmp);
               break;

            case (INTERVENTION_CLASS_RASTER):
               // Intervention class, is integer. If not an intervention, show INT_NODATA
               nTmp = static_cast<int>(pdScanline[nX]);

               if ((isnan(nTmp)) || (nTmp == m_nGISMissingValue))
               {
                  nTmp = m_nMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(nTmp);
               break;

            case (INTERVENTION_HEIGHT_RASTER):
               // Intervention height
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].SetInterventionHeight(dTmp);
               break;

            case (SUSP_SED_RASTER):
               // Initial Suspended Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].SetSuspendedSediment(dTmp);
               break;

            case (FINE_UNCONS_RASTER):
               // Initial Unconsolidated Fine Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dTmp);
               break;

            case (SAND_UNCONS_RASTER):
               // Initial Unconsolidated Sand Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dTmp);
               break;

            case (COARSE_UNCONS_RASTER):
               // Initial Unconsolidated Coarse Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dTmp);
               break;

            case (FINE_CONS_RASTER):
               // Initial Consolidated Fine Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->SetFineDepth(dTmp);
               break;

            case (SAND_CONS_RASTER):
               // Initial Consolidated Sand Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->SetSandDepth(dTmp);
               break;

            case (COARSE_CONS_RASTER):
               // Initial Consolidated Coarse Sediment GIS data
               dTmp = pdScanline[nX];

               if ((isnan(dTmp)) || (bFPIsEqual(dTmp, m_dGISMissingValue, TOLERANCE)))
               {
                  dTmp = m_dMissingValue;
                  nMissing++;
               }

               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->SetCoarseDepth(dTmp);
               break;
         }
      }
   }

   // Finished, so get rid of dataset object
   GDALClose(pGDALDataset);

   // Get rid of memory allocated to this array
   delete[] pdScanline;

   if (nMissing > 0)
   {
      cerr << WARN << nMissing << " missing values in " << strGISFile << endl;
      LogStream << WARN << nMissing << " missing values in " << strGISFile << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Writes GIS raster files using GDAL, using data from the RasterGrid array
//===============================================================================================================================
bool CSimulation::bWriteRasterGISFile(int const nDataItem, string const *strPlotTitle, int const nLayer, double const dElev)
{
   bool bIsInteger = false;
   bool bIsUnsignedLong = false;

   // Begin constructing the file name for this save
   string strFilePathName(m_strOutPath);
   string strLayer = "_layer_";

   stringstream ststrTmp;

   strLayer.append(to_string(nLayer + 1));

   switch (nDataItem)
   {
      case (RASTER_PLOT_BASEMENT_ELEVATION):
         strFilePathName.append(RASTER_BASEMENT_ELEVATION_NAME);
         break;

      case (RASTER_PLOT_SED_TOP_INC_TALUS_ELEV):
         strFilePathName.append(RASTER_SEDIMENT_TOP_ELEVATION_NAME);
         break;

      case (RASTER_PLOT_TOP_ELEV_INC_SEA):
         strFilePathName.append(RASTER_TOP_ELEVATION_INC_SEA_NAME);
         break;

      case (RASTER_PLOT_TALUS):
         strFilePathName.append(RASTER_TALUS_NAME);
         break;

      case (RASTER_PLOT_CONS_SED_SLOPE):
         strFilePathName.append(RASTER_SLOPE_OF_CONSOLIDATED_SEDIMENT_NAME);
         break;

      case (RASTER_PLOT_SLOPE_FOR_CLIFF_TOE):
         strFilePathName.append(RASTER_SLOPE_FOR_CLIFF_TOE_NAME);
         break;

      // case (RASTER_PLOT_CLIFF_TOE):
      //    strFilePathName.append(RASTER_CLIFF_TOE_NAME);
      //    break;

      case (RASTER_PLOT_SEA_DEPTH):
         strFilePathName.append(RASTER_SEA_DEPTH_NAME);
         break;

      case (RASTER_PLOT_AVG_SEA_DEPTH):
         strFilePathName.append(RASTER_AVG_SEA_DEPTH_NAME);
         break;

      case (RASTER_PLOT_WAVE_HEIGHT):
         strFilePathName.append(RASTER_WAVE_HEIGHT_NAME);
         break;

      case (RASTER_PLOT_AVG_WAVE_HEIGHT):
         strFilePathName.append(RASTER_AVG_WAVE_HEIGHT_NAME);
         break;

      case (RASTER_PLOT_WAVE_ORIENTATION):
         strFilePathName.append(RASTER_WAVE_ORIENTATION_NAME);
         break;

      case (RASTER_PLOT_AVG_WAVE_ORIENTATION):
         strFilePathName.append(RASTER_AVG_WAVE_ORIENTATION_NAME);
         break;

      case (RASTER_PLOT_BEACH_PROTECTION):
         strFilePathName.append(RASTER_BEACH_PROTECTION_NAME);
         break;

      case (RASTER_PLOT_POTENTIAL_PLATFORM_EROSION):
         strFilePathName.append(RASTER_POTENTIAL_PLATFORM_EROSION_NAME);
         break;

      case (RASTER_PLOT_ACTUAL_PLATFORM_EROSION):
         strFilePathName.append(RASTER_ACTUAL_PLATFORM_EROSION_NAME);
         break;

      case (RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION):
         strFilePathName.append(RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_NAME);
         break;

      case (RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION):
         strFilePathName.append(RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_NAME);
         break;

      case (RASTER_PLOT_POTENTIAL_BEACH_EROSION):
         strFilePathName.append(RASTER_POTENTIAL_BEACH_EROSION_NAME);
         break;

      case (RASTER_PLOT_ACTUAL_BEACH_EROSION):
         strFilePathName.append(RASTER_ACTUAL_BEACH_EROSION_NAME);
         break;

      case (RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION):
         strFilePathName.append(RASTER_TOTAL_POTENTIAL_BEACH_EROSION_NAME);
         break;

      case (RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION):
         strFilePathName.append(RASTER_TOTAL_ACTUAL_BEACH_EROSION_NAME);
         break;

      case (RASTER_PLOT_BEACH_DEPOSITION):
         strFilePathName.append(RASTER_BEACH_DEPOSITION_NAME);
         break;

      case (RASTER_PLOT_TOTAL_BEACH_DEPOSITION):
         strFilePathName.append(RASTER_TOTAL_BEACH_DEPOSITION_NAME);
         break;

      case (RASTER_PLOT_SUSPENDED_SEDIMENT):
         strFilePathName.append(RASTER_SUSP_SED_NAME);
         break;

      case (RASTER_PLOT_AVG_SUSPENDED_SEDIMENT):
         strFilePathName.append(RASTER_AVG_SUSP_SED_NAME);
         break;

      case (RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT):
         strFilePathName.append(RASTER_FINE_UNCONS_NAME);
         strFilePathName.append(strLayer);
         break;

      case (RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT):
         strFilePathName.append(RASTER_SAND_UNCONS_NAME);
         strFilePathName.append(strLayer);
         break;

      case (RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT):
         strFilePathName.append(RASTER_COARSE_UNCONS_NAME);
         strFilePathName.append(strLayer);
         break;

      case (RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT):
         strFilePathName.append(RASTER_FINE_CONS_NAME);
         strFilePathName.append(strLayer);
         break;

      case (RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT):
         strFilePathName.append(RASTER_SAND_CONS_NAME);
         strFilePathName.append(strLayer);
         break;

      case (RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT):
         strFilePathName.append(RASTER_COARSE_CONS_NAME);
         strFilePathName.append(strLayer);
         break;

      case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE):
         strFilePathName.append(RASTER_CLIFF_COLLAPSE_EROSION_FINE_NAME);
         break;

      case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND):
         strFilePathName.append(RASTER_CLIFF_COLLAPSE_EROSION_SAND_NAME);
         break;

      case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE):
         strFilePathName.append(RASTER_CLIFF_COLLAPSE_EROSION_COARSE_NAME);
         break;

      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE):
         strFilePathName.append(RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_NAME);
         break;

      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND):
         strFilePathName.append(RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_NAME);
         break;

      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE):
         strFilePathName.append(RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_NAME);
         break;

      case (RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND):
         strFilePathName.append(RASTER_CLIFF_COLLAPSE_DEPOSITION_SAND_NAME);
         break;

      case (RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE):
         strFilePathName.append(RASTER_CLIFF_COLLAPSE_DEPOSITION_COARSE_NAME);
         break;

      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND):
         strFilePathName.append(RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_NAME);
         break;

      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE):
         strFilePathName.append(RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_NAME);
         break;

      case (RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP):
         strFilePathName.append(RASTER_CLIFF_COLLAPSE_TIMESTEP_NAME);
         break;

      case (RASTER_PLOT_CLIFF_NOTCH_ALL):
         strFilePathName.append(RASTER_CLIFF_NOTCH_ALL_NAME);
         break;

      case (RASTER_PLOT_INTERVENTION_HEIGHT):
         strFilePathName.append(RASTER_INTERVENTION_HEIGHT_NAME);
         break;

      case (RASTER_PLOT_DEEP_WATER_WAVE_ORIENTATION):
         strFilePathName.append(RASTER_DEEP_WATER_WAVE_ORIENTATION_NAME);
         break;

      case (RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT):
         strFilePathName.append(RASTER_DEEP_WATER_WAVE_HEIGHT_NAME);
         break;

      case (RASTER_PLOT_POLYGON_GAIN_OR_LOSS):
         strFilePathName.append(RASTER_POLYGON_GAIN_OR_LOSS_NAME);
         break;

      case (RASTER_PLOT_DEEP_WATER_WAVE_PERIOD):
         strFilePathName.append(RASTER_WAVE_PERIOD_NAME);
         break;

      case (RASTER_PLOT_SEDIMENT_INPUT):
         strFilePathName.append(RASTER_SEDIMENT_INPUT_EVENT_NAME);
         break;

      case (RASTER_PLOT_BEACH_MASK):
         bIsInteger = true;
         strFilePathName.append(RASTER_BEACH_MASK_NAME);
         break;

      case (RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK):
         bIsInteger = true;
         strFilePathName.append(RASTER_POTENTIAL_PLATFORM_EROSION_MASK_NAME);
         break;

      case (RASTER_PLOT_INUNDATION_MASK):
         bIsInteger = true;
         strFilePathName.append(RASTER_INUNDATION_MASK_NAME);
         break;

      case (RASTER_PLOT_SLICE):
         bIsInteger = true;
         ststrTmp.str("");
         ststrTmp.clear();

         // TODO 031 Get working for multiple slices
         strFilePathName.append(RASTER_SLICE_NAME);
         ststrTmp << "_" << dElev << "_";
         strFilePathName.append(ststrTmp.str());
         break;

      case (RASTER_PLOT_LANDFORM):
         bIsInteger = true;
         strFilePathName.append(RASTER_LANDFORM_NAME);
         break;

      case (RASTER_PLOT_INTERVENTION_CLASS):
         bIsInteger = true;
         strFilePathName.append(RASTER_INTERVENTION_CLASS_NAME);
         break;

      case (RASTER_PLOT_COAST):
         bIsInteger = true;
         strFilePathName.append(RASTER_COAST_NAME);
         break;

      case (RASTER_PLOT_NORMAL_PROFILE):
         bIsInteger = true;
         strFilePathName.append(RASTER_COAST_NORMAL_NAME);
         break;

      case (RASTER_PLOT_ACTIVE_ZONE):
         bIsInteger = true;
         strFilePathName.append(RASTER_ACTIVE_ZONE_NAME);
         break;

      case (RASTER_PLOT_POLYGON):
         bIsInteger = true;
         strFilePathName.append(RASTER_POLYGON_NAME);
         break;

      case (RASTER_PLOT_SHADOW_ZONE):
         bIsInteger = true;
         strFilePathName.append(RASTER_SHADOW_ZONE_NAME);
         break;

      case (RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE):
         bIsInteger = true;
         strFilePathName.append(RASTER_SHADOW_DOWNDRIFT_ZONE_NAME);
         break;

      case (RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT):
         bIsInteger = true;
         strFilePathName.append(RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_NAME);
         break;

      // case (RASTER_PLOT_SETUP_SURGE_FLOOD_MASK):
      //    bIsInteger = true;
      //    strFilePathName.append(RASTER_SETUP_SURGE_FLOOD_MASK_NAME);
      //    break;
      //
      // case (RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK):
      //    bIsInteger = true;
      //    strFilePathName.append(RASTER_SETUP_SURGE_RUNUP_FLOOD_MASK_NAME);
      //    break;
      //
      // case (RASTER_PLOT_WAVE_FLOOD_LINE):
      //    bIsInteger = true;
      //    strFilePathName.append(RASTER_WAVE_FLOOD_LINE_NAME);
      //    break;
   }

   // Append the 'save number' to the filename, and prepend zeros to the save
   // number
   ststrTmp.str("");
   ststrTmp.clear();

   strFilePathName.append("_");

   if (m_bGISSaveDigitsSequential)
   {
      // Save number is m_bGISSaveDigitsSequential
      ststrTmp << FillToWidth('0', m_nGISMaxSaveDigits) << m_nGISSave;
   }
   else
   {
      // Save number is iteration
      ststrTmp << FillToWidth('0', m_nGISMaxSaveDigits) << m_ulIter;
   }

   strFilePathName.append(ststrTmp.str());

   // Finally, maybe append the extension
   if (! m_strGDALRasterOutputDriverExtension.empty())
   {
      strFilePathName.append(".");
      strFilePathName.append(m_strGDALRasterOutputDriverExtension);
   }

   // TODO 065 Used to try to debug floating point exception in pDriver->Create()
   // below CPLSetConfigOption("CPL_DEBUG", "ON");
   // CPLSetConfigOption("GDAL_NUM_THREADS", "1");

   GDALDriver *pDriver;
   GDALDataset *pDataSet;

   if (m_bGDALCanCreate)
   {
      // The user-requested raster driver supports the Create() method
      pDriver = GetGDALDriverManager()->GetDriverByName(
         m_strRasterGISOutFormat.c_str());

      if (bIsInteger)
      {
         pDataSet =
            pDriver->Create(strFilePathName.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Int16, m_papszGDALRasterOptions);
      }
      else if (bIsUnsignedLong)
      {
         pDataSet =
            pDriver->Create(strFilePathName.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_UInt32, m_papszGDALRasterOptions);
      }
      else if (m_strRasterGISOutFormat == "gpkg")
      {
         // TODO 065 Floating point exception here
         pDataSet = pDriver->Create(strFilePathName.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Byte, m_papszGDALRasterOptions);
      }
      else
      {
         pDataSet = pDriver->Create(strFilePathName.c_str(), m_nXGridSize, m_nYGridSize, 1, m_GDALWriteFloatDataType, m_papszGDALRasterOptions);
      }

      if (NULL == pDataSet)
      {
         // Error, couldn't create file
         cerr << ERR << "cannot create " << m_strRasterGISOutFormat << " file named " << strFilePathName << endl;
         return false;
      }
   }
   else
   {
      // The user-requested raster driver does not support the Create() method, so we must first create a memory-file dataset
      pDriver = GetGDALDriverManager()->GetDriverByName("MEM");
      pDataSet = pDriver->Create("", m_nXGridSize, m_nYGridSize, 1, m_GDALWriteFloatDataType, NULL);

      if (NULL == pDataSet)
      {
         // Couldn't create in-memory file dataset
         cerr << ERR << "cannot create in-memory file for " << m_strRasterGISOutFormat << " file named " << strFilePathName << endl << CPLGetLastErrorMsg() << endl;
         return false;
      }
   }

   // Set projection info for output dataset (will be same as was read in from basement DEM)
   CPLPushErrorHandler(CPLQuietErrorHandler);                            // Needed to get next line to fail silently, if it fails
   pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());      // Will fail for some formats
   CPLPopErrorHandler();

   // Set geotransformation info for output dataset (will be same as was read in from DEM)
   if (CE_Failure == pDataSet->SetGeoTransform(m_dGeoTransform))
      LogStream << WARN << "cannot write geotransformation information to " << m_strRasterGISOutFormat << " file named " << strFilePathName << endl << CPLGetLastErrorMsg() << endl;

   // Allocate memory for a 1D array, to hold the floating point raster band data for GDAL
   double *pdRaster = new double[m_ulNumCells];
   if (NULL == pdRaster)
   {
      // Error, can't allocate memory
      cerr << ERR << "cannot allocate memory for " << m_ulNumCells << " x 1D floating-point array for " << m_strRasterGISOutFormat << " file named " << strFilePathName << endl;
      return (RTN_ERR_MEMALLOC);
   }

   bool bScaleOutput = false;
   double dRangeScale = 0;
   double dDataMin = 0;

   if (! m_bGDALCanWriteFloat)
   {
      double dDataMax = 0;

      // The output file format cannot handle floating-point numbers, so we may need to scale the output
      GetRasterOutputMinMax(nDataItem, dDataMin, dDataMax, nLayer, 0);

      double const dDataRange = dDataMax - dDataMin;
      double const dWriteRange = static_cast<double>(m_lGDALMaxCanWrite - m_lGDALMinCanWrite);

      if (dDataRange > 0)
         dRangeScale = dWriteRange / dDataRange;

      // If we are attempting to write values which are outside this format's allowable range, and the user has set the option, then scale the output
      if (((dDataMin < static_cast<double>(m_lGDALMinCanWrite)) || (dDataMax > static_cast<double>(m_lGDALMaxCanWrite))) && m_bScaleRasterOutput)
         bScaleOutput = true;
   }

   // Fill the array
   int n = 0;
   int nPoly = 0;
   int nPolyCoast = 0;
   int nTopLayer = 0;
   double dTmp = 0;

   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         switch (nDataItem)
         {
            case (RASTER_PLOT_BASEMENT_ELEVATION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBasementElev();
               break;

            case (RASTER_PLOT_SED_TOP_INC_TALUS_ELEV):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus();
               break;

            case (RASTER_PLOT_TOP_ELEV_INC_SEA):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTopElevIncSea();
               break;

            case (RASTER_PLOT_TALUS):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTalusDepth();
               break;

            case (RASTER_PLOT_CONS_SED_SLOPE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedSlope();
               break;

            // case (RASTER_PLOT_SLOPE_FOR_CLIFF_TOE):
            //    dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSlopeForCliffToe();
            //    break;
            //
            // case (RASTER_PLOT_CLIFF_TOE):
            //    dTmp = static_cast<double>(m_pRasterGrid->m_Cell[nX][nY].bIsCliffToe());
            //    break;

            case (RASTER_PLOT_SEA_DEPTH):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();
               break;

            case (RASTER_PLOT_AVG_SEA_DEPTH):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSeaDepth() / static_cast<double>(m_ulIter);
               break;

            case (RASTER_PLOT_WAVE_HEIGHT):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_AVG_WAVE_HEIGHT):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveHeight() / static_cast<double>(m_ulIter);
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_WAVE_ORIENTATION):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_AVG_WAVE_ORIENTATION):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveAngle() / static_cast<double>(m_ulIter);
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_BEACH_PROTECTION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor();

               if (bFPIsEqual(dTmp, DBL_NODATA, TOLERANCE))
                  dTmp = m_dMissingValue;
               else
                  dTmp = 1 - dTmp;      // Output the inverse, seems more intuitive
               break;

            case (RASTER_PLOT_POTENTIAL_PLATFORM_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion();
               break;

            case (RASTER_PLOT_ACTUAL_PLATFORM_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetActualPlatformErosion();
               break;

            case (RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotPotentialPlatformErosion();
               break;

            case (RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotActualPlatformErosion();
               break;

            case (RASTER_PLOT_POTENTIAL_BEACH_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetPotentialBeachErosion();
               break;

            case (RASTER_PLOT_ACTUAL_BEACH_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetActualBeachErosion();
               break;

            case (RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotPotentialBeachErosion();
               break;

            case (RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotActualBeachErosion();
               break;

            case (RASTER_PLOT_BEACH_DEPOSITION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBeachDeposition();
               break;

            case (RASTER_PLOT_TOTAL_BEACH_DEPOSITION):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotBeachDeposition();
               break;

            case (RASTER_PLOT_SUSPENDED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();
               break;

            case (RASTER_PLOT_AVG_SUSPENDED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSuspendedSediment() / static_cast<double>(m_ulIter);
               break;

            case (RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
               break;

            case (RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
               break;

            case (RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
               break;

            case (RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->dGetFineDepth();
               break;

            case (RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->dGetSandDepth();
               break;

            case (RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();
               break;

            case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionFine();
               break;

            case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionSand();
               break;

            case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionCoarse();
               break;

            case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapseFine();
               break;

            case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapseSand();
               break;

            case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapseCoarse();
               break;

            case (RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseSandTalusDeposition();
               break;

            case (RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseCoarseTalusDeposition();
               break;

            case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSandTalusDeposition();
               break;

            case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotCoarseTalusDeposition();
               break;

            case (RASTER_PLOT_CLIFF_NOTCH_ALL):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetCliffNotchIncisionDepth();
               break;

            case (RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP):
               dTmp = static_cast<double>(m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->ulGetCliffCollapseTimestep());
               bIsUnsignedLong = true;
               break;

            case (RASTER_PLOT_INTERVENTION_HEIGHT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight();
               break;

            case (RASTER_PLOT_DEEP_WATER_WAVE_ORIENTATION):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle();
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight();
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_DEEP_WATER_WAVE_PERIOD):
               if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
                  dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWavePeriod();
               else
                  dTmp = 0;
               break;

            case (RASTER_PLOT_POLYGON_GAIN_OR_LOSS):
               nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
               nPolyCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();

               if (nPoly == INT_NODATA)
                  dTmp = m_dMissingValue;
               else
               {
                  // Get total volume (all sediment size classes) of change in sediment for this polygon for this timestep (-ve erosion, +ve deposition)
                  dTmp = m_VCoast[nPolyCoast].pGetPolygon(nPoly)->dGetBeachDepositionAndSuspensionAllUncons() * m_dCellArea;

                  // Calculate the rate in m^3 / sec
                  dTmp /= (m_dTimeStep * 3600);
               }
               break;

            case (RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK):
               // cppcheck-suppress assignBoolToFloat
               dTmp = m_pRasterGrid->m_Cell[nX][nY].bPotentialPlatformErosion();
               break;

            case (RASTER_PLOT_INUNDATION_MASK):
               // cppcheck-suppress assignBoolToFloat
               dTmp = m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea();
               break;

            case (RASTER_PLOT_BEACH_MASK):
               dTmp = 0;
               nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

               if ((nTopLayer == INT_NODATA) || (nTopLayer == NO_NONZERO_THICKNESS_LAYERS))
                  break;

               if ((m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->dGetAllUnconsDepth() > 0) && (m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() > m_dThisIterSWL))
                  dTmp = 1;

               break;

            case (RASTER_PLOT_SLICE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dElev);
               break;

            case (RASTER_PLOT_LANDFORM):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory();
               bIsInteger = true;
               break;

            case (RASTER_PLOT_INTERVENTION_CLASS):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetInterventionClass();
               bIsInteger = true;
               break;

            case (RASTER_PLOT_COAST):
               dTmp = (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() ? 1 : 0);
               break;

            case (RASTER_PLOT_NORMAL_PROFILE):
               // dTmp = (m_pRasterGrid->m_Cell[nX][nY].bIsProfile() ? 1 : 0);
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();
               bIsInteger = true;
               break;

            case (RASTER_PLOT_ACTIVE_ZONE):
               dTmp = (m_pRasterGrid->m_Cell[nX][nY].bIsInActiveZone() ? 1 : 0);
               break;

            case (RASTER_PLOT_POLYGON):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
               bIsInteger = true;
               break;

            case (RASTER_PLOT_SHADOW_ZONE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();
               bIsInteger = true;
               break;

            case (RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
               bIsInteger = true;
               break;

            case (RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT):
               nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
               nPolyCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();
               bIsInteger = true;

               if (nPoly == INT_NODATA)
                  dTmp = m_nMissingValue;
               else
               {
                  if (m_VCoast[nPolyCoast].pGetPolygon(nPoly)->bDownCoastThisIter())
                     dTmp = 1;
                  else
                     dTmp = 0;
               }
               break;

            case (RASTER_PLOT_SEDIMENT_INPUT):
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetTotAllSedimentInputDepth();
               break;

            // case (RASTER_PLOT_SETUP_SURGE_FLOOD_MASK):
            //    dTmp = (m_pRasterGrid->m_Cell[nX][nY].bIsFloodBySetupSurge() ? 1 : 0);
            //    break;
            //
            // case (RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK):
            //    dTmp = (m_pRasterGrid->m_Cell[nX][nY].bIsFloodBySetupSurgeRunup() ? 1 : 0);
            //    break;

            // case (RASTER_PLOT_WAVE_FLOOD_LINE):
            //    dTmp = (m_pRasterGrid->m_Cell[nX][nY].bIsFloodline() ? 1 : 0);
            //    break;
         }

         // If necessary, scale this value
         if (bScaleOutput)
         {
            if (bFPIsEqual(dTmp, DBL_NODATA, TOLERANCE))
               dTmp = 0;      // TODO 032 Improve this
            else
               dTmp = dRound(static_cast<double>(m_lGDALMinCanWrite) + (dRangeScale * (dTmp - dDataMin)));
         }

         // Write this value to the array
         pdRaster[n++] = dTmp;
      }
   }

   // Create a single raster band
   GDALRasterBand *pBand = pDataSet->GetRasterBand(1);

   // And fill it with the NODATA value
   if (bIsInteger)
      pBand->Fill(m_nMissingValue);
   else if (bIsUnsignedLong)
      pBand->Fill(static_cast<double>(m_ulMissingValue));
   else
      pBand->Fill(m_dMissingValue);

   // Set value units for this band
   string strUnits;

   switch (nDataItem)
   {
      case (RASTER_PLOT_ACTUAL_BEACH_EROSION):
      case (RASTER_PLOT_ACTUAL_PLATFORM_EROSION):
      case (RASTER_PLOT_AVG_SEA_DEPTH):
      case (RASTER_PLOT_AVG_SUSPENDED_SEDIMENT):
      case (RASTER_PLOT_AVG_WAVE_HEIGHT):
      case (RASTER_PLOT_BASEMENT_ELEVATION):
      case (RASTER_PLOT_BEACH_DEPOSITION):
      case (RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE):
      case (RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND):
      case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE):
      case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE):
      case (RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND):
      case (RASTER_PLOT_CLIFF_NOTCH_ALL):
      case (RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT):
      case (RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT):
      case (RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT):
      case (RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT):
      case (RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT):
      case (RASTER_PLOT_INTERVENTION_HEIGHT):
      case (RASTER_PLOT_POTENTIAL_BEACH_EROSION):
      case (RASTER_PLOT_POTENTIAL_PLATFORM_EROSION):
      case (RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT):
      case (RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT):
      case (RASTER_PLOT_SEA_DEPTH):
      case (RASTER_PLOT_SEDIMENT_INPUT):
      case (RASTER_PLOT_SED_TOP_INC_TALUS_ELEV):
      case (RASTER_PLOT_SUSPENDED_SEDIMENT):
      case (RASTER_PLOT_TALUS):
      case (RASTER_PLOT_TOP_ELEV_INC_SEA):
      case (RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION):
      case (RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION):
      case (RASTER_PLOT_TOTAL_BEACH_DEPOSITION):
      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE):
      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND):
      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE):
      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE):
      case (RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND):
      case (RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION):
      case (RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION):
      case (RASTER_PLOT_WAVE_HEIGHT):
         strUnits = "m";
         break;

      case (RASTER_PLOT_CONS_SED_SLOPE):
         strUnits = "m/m";
         break;

         strUnits = "degrees";
      case (RASTER_PLOT_AVG_WAVE_ORIENTATION):
      case (RASTER_PLOT_WAVE_ORIENTATION):
         break;

      case (RASTER_PLOT_POLYGON_GAIN_OR_LOSS):
         strUnits = "cumecs";
         break;

      case (RASTER_PLOT_DEEP_WATER_WAVE_PERIOD):
         strUnits = "secs";
         break;

         strUnits = "none";
      case (RASTER_PLOT_ACTIVE_ZONE):
      case (RASTER_PLOT_BEACH_MASK):
      case (RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP):
      case (RASTER_PLOT_COAST):
      case (RASTER_PLOT_INTERVENTION_CLASS):
      case (RASTER_PLOT_INUNDATION_MASK):
      case (RASTER_PLOT_LANDFORM):
      case (RASTER_PLOT_NORMAL_PROFILE):
      case (RASTER_PLOT_POLYGON):
      case (RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT):
      case (RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK):
      // case (RASTER_PLOT_SETUP_SURGE_FLOOD_MASK):
      // case (RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK):
      case (RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE):
      case (RASTER_PLOT_SHADOW_ZONE):
      case (RASTER_PLOT_SLICE):
      case (RASTER_PLOT_WAVE_FLOOD_LINE):
         break;
   }

   CPLPushErrorHandler(CPLQuietErrorHandler);            // Needed to get next line to fail silently, if it fails
   pBand->SetUnitType(strUnits.c_str());                 // Not supported for some GIS formats
   CPLPopErrorHandler();

   // Tell the output dataset about NODATA (missing values)
   CPLPushErrorHandler(CPLQuietErrorHandler);            // Needed to get next line to fail silently, if it fails

   if (bIsInteger)
      pBand->SetNoDataValue(m_nMissingValue);            // Will fail for some formats
   if (bIsUnsignedLong)
      pBand->SetNoDataValueAsUInt64(m_ulMissingValue);   // Will fail for some formats
   else
      pBand->SetNoDataValue(m_dMissingValue);            // Will fail for some formats

   CPLPopErrorHandler();

   // Construct the description
   string strDesc(*strPlotTitle);

   if (nDataItem == RASTER_PLOT_SLICE)
   {
      ststrTmp.clear();
      ststrTmp << dElev << "m, ";
      strDesc.append(ststrTmp.str());
   }

   strDesc.append(" at ");
   strDesc.append(strDispTime(m_dSimElapsed, false, false));

   // Set the GDAL description
   pBand->SetDescription(strDesc.c_str());

   // Set raster category names
   char **papszCategoryNames = NULL;

   switch (nDataItem)
   {
      case (RASTER_PLOT_SLICE):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Basement");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 0");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 1");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 2");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 3");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 4");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 5");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 6");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 7");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 8");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Layer 9");
         break;

      case (RASTER_PLOT_LANDFORM):
         papszCategoryNames = CSLAddString(papszCategoryNames, "None");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Hinterland");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Sea");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Cliff");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Drift");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Intervention");

         papszCategoryNames = CSLAddString(papszCategoryNames, "Cliff on Coastline");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Inland Cliff");

         papszCategoryNames = CSLAddString(papszCategoryNames, "Mixed Drift");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Talus");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Beach");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Dunes");
         break;

      case (RASTER_PLOT_INTERVENTION_CLASS):
         papszCategoryNames = CSLAddString(papszCategoryNames, "None");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Structural");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Non-Structural");
         break;

      case (RASTER_PLOT_COAST):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Not coastline");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Coastline");
         break;

      case (RASTER_PLOT_NORMAL_PROFILE):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Not coastline-normal profile");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Coastline-normal profile");
         break;

      case (RASTER_PLOT_ACTIVE_ZONE):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Not in active zone");
         papszCategoryNames = CSLAddString(papszCategoryNames, "In active zone");
         break;

      case (RASTER_PLOT_POLYGON):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Not polygon");
         papszCategoryNames = CSLAddString(papszCategoryNames, "In polygon");
         break;

      case (RASTER_PLOT_SHADOW_ZONE):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Not in shadow zone");
         papszCategoryNames = CSLAddString(papszCategoryNames, "In shadow zone");
         break;

      case (RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Not in shadow downdrift zone");
         papszCategoryNames = CSLAddString(papszCategoryNames, "In shadow downdrift zone");
         break;

      case (RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT):
         papszCategoryNames = CSLAddString(papszCategoryNames, "Updrift movement of unconsolidated sediment ");
         papszCategoryNames = CSLAddString(papszCategoryNames, "Downdrift movement of unconsolidated sediment");
         break;

      // case (RASTER_PLOT_SETUP_SURGE_FLOOD_MASK):
      //    papszCategoryNames = CSLAddString(papszCategoryNames, "Inundated by swl setup and surge ");
      //    papszCategoryNames = CSLAddString(papszCategoryNames, "Not inundated by swl setup and surge");
      //    break;
      //
      // case (RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK):
      //    papszCategoryNames = CSLAddString(papszCategoryNames, "Inundated by swl setup, surge and runup ");
      //    papszCategoryNames = CSLAddString(papszCategoryNames, "Not inundated by swl setup, surge and runup");
      //    break;
      //
      // case (RASTER_PLOT_WAVE_FLOOD_LINE):
      //    papszCategoryNames = CSLAddString(papszCategoryNames, "Intersection line of inundation ");
      //    papszCategoryNames = CSLAddString(papszCategoryNames, "Not inundated by swl waves and runup");
      //    break;
   }

   CPLPushErrorHandler(CPLQuietErrorHandler);        // Needed to get next line to fail silently, if it fails
   pBand->SetCategoryNames(papszCategoryNames);      // Not supported for some GIS formats
   CPLPopErrorHandler();

   // Now write the data with optimized I/O
   // Enable multi-threaded compression for faster writing
   CPLSetThreadLocalConfigOption("GDAL_NUM_THREADS", "ALL_CPUS");

   if (CE_Failure == pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL))
   {
      // Write error, better error message
      cerr << ERR << "cannot write data for " << m_strRasterGISOutFormat << " file named " << strFilePathName << endl << CPLGetLastErrorMsg() << endl;
      delete[] pdRaster;
      return false;
   }

   // Calculate statistics for this band
   double dMin, dMax, dMean, dStdDev;
   CPLPushErrorHandler(CPLQuietErrorHandler);      // Needed to get next line to fail silently, if it fails
   pBand->ComputeStatistics(false, &dMin, &dMax, &dMean, &dStdDev, NULL, NULL);
   CPLPopErrorHandler();

   // And then write the statistics
   CPLPushErrorHandler(CPLQuietErrorHandler);      // Needed to get next line to fail silently, if it fails
   pBand->SetStatistics(dMin, dMax, dMean, dStdDev);
   CPLPopErrorHandler();

   if (! m_bGDALCanCreate)
   {
      // Since the user-selected raster driver cannot use the Create() method, we have been writing to a dataset created by the in-memory driver. So now we need to use CreateCopy() to copy this in-memory dataset to a file in the user-specified raster driver format
      GDALDriver *pOutDriver = GetGDALDriverManager()->GetDriverByName(m_strRasterGISOutFormat.c_str());
      GDALDataset *pOutDataSet = pOutDriver->CreateCopy(strFilePathName.c_str(), pDataSet, false, m_papszGDALRasterOptions, NULL, NULL);

      if (NULL == pOutDataSet)
      {
         // Couldn't create file
         cerr << ERR << "cannot create " << m_strRasterGISOutFormat << " file named " << strFilePathName << endl << CPLGetLastErrorMsg() << endl;
         return false;
      }

      // Get rid of this user-selected dataset object
      GDALClose(pOutDataSet);
   }

   // Get rid of dataset object
   GDALClose(pDataSet);

   // Also get rid of memory allocated to this array
   delete[] pdRaster;

   return true;
}

//===============================================================================================================================
//! Interpolates wave properties from profile points to all cells within polygons
//!
//! ALGORITHM: k-Nearest Neighbor Inverse Distance Weighting (k-NN IDW)
//!
//! This function takes wave height components (X and Y) measured at discrete profile points
//! and interpolates them to a regular grid using a spatial interpolation method.
//!
//! METHOD OVERVIEW:
//! 1. Builds a k-d tree spatial index from input profile point coordinates (X, Y)
//! 2. For each grid cell, finds the k nearest profile points (typically k=12)
//! 3. Calculates interpolated value using inverse distance weighting (IDW) with power=2
//! 4. Converts X/Y wave components back to magnitude and direction
//! 5. Updates grid cells with interpolated wave properties
//!
//! KEY PARAMETERS (defined in spatial_interpolation.cpp DualSpatialInterpolator):
//! - k_neighbors = 12    Number of nearest points to use for interpolation
//! - power = 2.0         Exponent for inverse distance weighting
//!                       (higher = more weight to closer points)
//!
//! TUNING GUIDANCE:
//! - Increase k_neighbors (e.g., 15-20) for smoother results with more averaging
//! - Decrease k_neighbors (e.g., 8-10) for results that follow local variations more closely
//! - Increase power (e.g., 3.0-4.0) to emphasize nearby points (sharper transitions)
//! - Decrease power (e.g., 1.0-1.5) for smoother, more gradual transitions
//!
//! IMPORTANT NOTES:
//! - This method does NOT respect transect structure (point 1 vs point 2, etc.)
//! - Treats all input points equally regardless of which transect they belong to
//! - Works well for scattered points but may not preserve transect-aligned features
//! - For transect-aware interpolation, consider bilinear methods instead
//!
//! COASTLINE ORIENTATION:
//! - The coastlines are available in m_VCoast[] vector
//! - Each coast has flux orientation: m_VCoast[i].dGetFluxOrientation(coastpoint)
//! - Each coast has breaking wave angle: m_VCoast[i].dGetBreakingWaveAngle(coastpoint)
//! - These could be used to implement coast-aware interpolation algorithms
//!
//! @param pVTransects         Vector of TransectWaveData containing wave data per transect
//! @param pVdDeepWaterX       X coordinates of deep water grid edge points
//! @param pVdDeepWaterY       Y coordinates of deep water grid edge points
//! @param pVdDeepWaterHeightX X component of wave height at deep water points
//! @param pVdDeepWaterHeightY Y component of wave height at deep water points
//! @return RTN_OK on success, error code otherwise
//===============================================================================================================================
int CSimulation::nInterpolateWavesToPolygonCells(vector<TransectWaveData> const *pVTransects, vector<double> const *pVdDeepWaterX, vector<double> const *pVdDeepWaterY, vector<double> const *pVdDeepWaterHeightX,  vector<double> const *pVdDeepWaterHeightY)
{
   // ============================================================================
   // STEP 1: Calculate grid dimensions and initialize variables
   // ============================================================================

   int nXSize = 0;
   int nYSize = 0;

   // Average values used as fallback when interpolation fails or returns NaN
   double dXAvg = 0;
   double dYAvg = 0;

   // Calculate bounding box size
   nXSize = m_nXMaxBoundingBox - m_nXMinBoundingBox + 1;
   nYSize = m_nYMaxBoundingBox - m_nYMinBoundingBox + 1;
   int const nGridSize = nXSize * nYSize;

   // Count total points across all transects plus deep water points
   unsigned int nPoints = 0;
   for (auto const &transect : *pVTransects)
   {
      nPoints += static_cast<unsigned int>(transect.VdX.size());
   }
   nPoints += static_cast<unsigned int>(pVdDeepWaterX->size());

   // Initialize output arrays (will hold interpolated X and Y wave components)
   vector<double> VdOutX(nGridSize, 0);
   vector<double> VdOutY(nGridSize, 0);

   // ============================================================================
   // STEP 2: Prepare input data for spatial interpolation
   // ============================================================================

   // Flatten transect data and deep water data into contiguous arrays for the interpolator
   std::vector<Point2D> points;
   std::vector<double> VdHeightX;
   std::vector<double> VdHeightY;

   points.reserve(nPoints);
   VdHeightX.reserve(nPoints);
   VdHeightY.reserve(nPoints);

   // Add profile/transect points
   for (auto const &transect : *pVTransects)
   {
      for (size_t i = 0; i < transect.VdX.size(); i++)
      {
         points.emplace_back(transect.VdX[i], transect.VdY[i]);
         VdHeightX.push_back(transect.VdHeightX[i]);
         VdHeightY.push_back(transect.VdHeightY[i]);
      }
   }

   // Add deep water grid edge points
   for (size_t i = 0; i < pVdDeepWaterX->size(); i++)
   {
      points.emplace_back((*pVdDeepWaterX)[i], (*pVdDeepWaterY)[i]);
      VdHeightX.push_back((*pVdDeepWaterHeightX)[i]);
      VdHeightY.push_back((*pVdDeepWaterHeightY)[i]);
   }

   // ============================================================================
   // STEP 3: Create spatial interpolator
   // ============================================================================
   //
   // DualSpatialInterpolator parameters:
   // - points: Input point coordinates (from profiles/transects)
   // - VdHeightX, VdHeightY: Wave height X and Y components at those points
   // - k_neighbors = 12: Use 12 nearest neighbors for interpolation
   //                     ** ADJUST THIS to change smoothness vs local detail **
   // - power = 2.0: Inverse distance weighting power
   //                ** ADJUST THIS to change influence of nearby vs distant points **
   //
   // The interpolator builds a k-d tree for fast nearest neighbor search and shares it between X and Y interpolation for efficiency
   DualSpatialInterpolator interp(points, VdHeightX, VdHeightY, 6, 2.0);

   // ============================================================================
   // STEP 4: Build query points (grid cells where we want interpolated values)
   // ============================================================================

   std::vector<Point2D> query_points;
   query_points.reserve(nGridSize);
   for (int nY = m_nYMinBoundingBox; nY <= m_nYMaxBoundingBox; nY++)
   {
      for (int nX = m_nXMinBoundingBox; nX <= m_nXMaxBoundingBox; nX++)
      {
         query_points.emplace_back(static_cast<double>(nX), static_cast<double>(nY));
      }
   }

   // ============================================================================
   // STEP 5: Perform batch interpolation
   // ============================================================================
   //
   // This does the actual interpolation for all grid points at once. Uses OpenMP parallelization if available (see spatial_interpolation.cpp). Interpolates both X and Y components simultaneously using shared k-d tree
   interp.Interpolate(query_points, VdOutX, VdOutY);

   // ============================================================================
   // STEP 6: Validate results and calculate average values for fallback
   // ============================================================================
   //
   // Check for NaN or unreasonably large values and replace with missing value marker
   // Also calculate average of valid values to use as fallback

   int nXValid = 0;
   int nYValid = 0;

   // Validate X component
   for (unsigned int n = 0; n < VdOutX.size(); n++)
   {
      if (isnan(VdOutX[n]))
         VdOutX[n] = m_dMissingValue;
      else if (tAbs(VdOutX[n]) > 1e10)      // Sanity check for unreasonably large values
         VdOutX[n] = m_dMissingValue;
      else
      {
         dXAvg += VdOutX[n];
         nXValid++;
      }
   }

   // Validate Y component
   for (unsigned int n = 0; n < VdOutY.size(); n++)
   {
      if (isnan(VdOutY[n]))
         VdOutY[n] = m_dMissingValue;
      else if (tAbs(VdOutY[n]) > 1e10)      // Sanity check for unreasonably large values
         VdOutY[n] = m_dMissingValue;
      else
      {
         dYAvg += VdOutY[n];
         nYValid++;
      }
   }

   // Calculate averages (for use as fallback when individual cells have missing values)
   if (nXValid > 0)
      dXAvg /= nXValid;
   if (nYValid > 0)
      dYAvg /= nYValid;

   // ============================================================================
   // STEP 7: Update grid cells with interpolated wave properties
   // ============================================================================
   //
   // Convert X and Y components back to magnitude and direction, then update each cell's wave attributes

   int n = 0;

   for (int nY = 0; nY < nYSize; nY++)
   {
      for (int nX = 0; nX < nXSize; nX++)
      {
         int const nActualX = nX + m_nXMinBoundingBox;
         int const nActualY = nY + m_nYMinBoundingBox;

         if (m_pRasterGrid->m_Cell[nActualX][nActualY].bIsInContiguousSea())
         {
            // Only update sea cells

            if (m_pRasterGrid->m_Cell[nActualX][nActualY].nGetPolygonID() == INT_NODATA)
            {
               // --------------------------------------------------------------------
               // Deep water cell (NOT in a polygon)
               // --------------------------------------------------------------------
               // Use the cell's pre-assigned deep water wave values (these cells are beyond the coastal zone, so don't need interpolation)

               double const dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[nActualX][nActualY].dGetCellDeepWaterWaveHeight();
               m_pRasterGrid->m_Cell[nActualX][nActualY].SetWaveHeight(dDeepWaterWaveHeight);

               double const dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[nActualX][nActualY].dGetCellDeepWaterWaveAngle();
               m_pRasterGrid->m_Cell[nActualX][nActualY].SetWaveAngle(dDeepWaterWaveAngle);
            }
            else
            {
               // --------------------------------------------------------------------
               // Coastal zone cell (IN a polygon)
               // --------------------------------------------------------------------
               // Use the interpolated wave values calculated above

               double dWaveHeightX;
               double dWaveHeightY;

               // Get interpolated X component (use average as fallback if missing/invalid)
               if ((isnan(VdOutX[n])) || (bFPIsEqual(VdOutX[n], m_dMissingValue, TOLERANCE)))
                  dWaveHeightX = dXAvg;
               else
                  dWaveHeightX = VdOutX[n];

               // Get interpolated Y component (use average as fallback if missing/invalid)
               if ((isnan(VdOutY[n])) || (bFPIsEqual(VdOutY[n], m_dMissingValue, TOLERANCE)))
                  dWaveHeightY = dYAvg;
               else
                  dWaveHeightY = VdOutY[n];

               // Convert X/Y components to magnitude and direction
               double const dWaveHeight = sqrt((dWaveHeightX * dWaveHeightX) + (dWaveHeightY * dWaveHeightY));
               double const dWaveDir = atan2(dWaveHeightX, dWaveHeightY) * (180 / PI);

               // Update the cell's wave attributes
               m_pRasterGrid->m_Cell[nActualX][nActualY].SetWaveHeight(
                  dWaveHeight);
               m_pRasterGrid->m_Cell[nActualX][nActualY].SetWaveAngle(
                  dKeepWithin360(dWaveDir));

               // Calculate wave height-to-depth ratio and update active zone status (active zone = where waves are breaking or near-breaking)
               double const dSeaDepth =
                  m_pRasterGrid->m_Cell[nActualX][nActualY].dGetSeaDepth();

               if ((dWaveHeight / dSeaDepth) >= m_dBreakingWaveHeightDepthRatio)
                  m_pRasterGrid->m_Cell[nActualX][nActualY].SetInActiveZone(true);

               // LogStream << " nX = " << nX << " nY = " << nY << " [" <<
               // nActualX
               // << "][" << nActualY << "] waveheight = " << dWaveHeight << "
               // dWaveDir = " << dWaveDir << " dKeepWithin360(dWaveDir) = " <<
               // dKeepWithin360(dWaveDir) << endl;
            }
         }

         // Increment with safety check
         n++;
         n = tMin(n, static_cast<int>(VdOutX.size() - 1));
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! If the user supplies multiple deep water wave height and angle values, this routine interplates these to all cells (including dry land cells)
//===============================================================================================================================
int CSimulation::nInterpolateAllDeepWaterWaveValues(void)
{
   // Interpolate deep water height and orientation from multiple user-supplied values
   unsigned int const nUserPoints = static_cast<unsigned int>(m_VdDeepWaterWaveStationX.size());

   // Performance optimization: Enable GDAL threading for interpolation
   CPLSetThreadLocalConfigOption("GDAL_NUM_THREADS", "ALL_CPUS");

   // Call GDALGridCreate() with the GGA_InverseDistanceToAPower interpolation algorithm. It has following parameters: radius1 is the first radius (X axis if rotation angle is 0) of the search ellipse, set this to zero (the default) to use the whole point array; radius2 is the second radius (Y axis if rotation angle is 0) of the search ellipse, again set this parameter to zero (the default) to use the whole point array; angle is the angle of the search ellipse rotation in degrees (counter clockwise, default 0.0); nodata is the NODATA marker to fill empty points (default 0.0) TODO 086
   GDALGridInverseDistanceToAPowerOptions *pOptions = new GDALGridInverseDistanceToAPowerOptions();
   pOptions->dfAngle = 0;
   pOptions->dfAnisotropyAngle = 0;
   pOptions->dfAnisotropyRatio = 0;
   pOptions->dfPower = 2;              // Reduced from 3 to 2 for faster computation
   pOptions->dfSmoothing = 50;         // Reduced from 100 to 50 for faster computation
   pOptions->dfRadius1 = 0;
   pOptions->dfRadius2 = 0;
   pOptions->nMaxPoints = 12;          // Limit points for faster computation (was 0 = unlimited)
   pOptions->nMinPoints = 3;           // Minimum points needed for interpolation
   pOptions->dfNoDataValue = m_nMissingValue;

   // CPLSetConfigOption("CPL_DEBUG", "ON");
   // CPLSetConfigOption("GDAL_NUM_THREADS", "1");

   // OK, now create a gridded version of wave height: first create the GDAL context TODO 086
   // GDALGridContext* pContext = GDALGridContextCreate(GGA_InverseDistanceToAPower, pOptions, nUserPoints, &m_VdDeepWaterWaveStationX[0], &m_VdDeepWaterWaveStationY[0], &m_VdThisIterDeepWaterWaveStationHeight[0], true);
   GDALGridContext *pContext = GDALGridContextCreate(GGA_InverseDistanceToAPower, pOptions, nUserPoints, m_VdDeepWaterWaveStationX.data(), m_VdDeepWaterWaveStationY.data(), m_VdThisIterDeepWaterWaveStationHeight.data(), true);

   if (pContext == NULL)
   {
      delete pOptions;
      return RTN_ERR_GRIDCREATE;
   }

   // Now process the context
   double *dHeightOut = new double[m_ulNumCells];
   int nRet = GDALGridContextProcess(pContext, 0, m_nXGridSize - 1, 0, m_nYGridSize - 1, m_nXGridSize, m_nYGridSize, GDT_Float64, dHeightOut, NULL, NULL);

   if (nRet == CE_Failure)
   {
      delete[] dHeightOut;
      delete pOptions;
      return RTN_ERR_GRIDCREATE;
   }

   // Get rid of the context
   GDALGridContextFree(pContext);

   // Next create a gridded version of wave orientation: first create the GDAL context
   // pContext = GDALGridContextCreate(GGA_InverseDistanceToAPower, pOptions, nUserPoints,  &(m_VdDeepWaterWaveStationX[0]), &(m_VdDeepWaterWaveStationY[0]), (&m_VdThisIterDeepWaterWaveStationAngle[0]), true);
   pContext = GDALGridContextCreate(GGA_InverseDistanceToAPower, pOptions, nUserPoints, m_VdDeepWaterWaveStationX.data(), m_VdDeepWaterWaveStationY.data(), m_VdThisIterDeepWaterWaveStationAngle.data(), true);

   if (pContext == NULL)
   {
      delete[] dHeightOut;
      delete pOptions;
      return RTN_ERR_GRIDCREATE;
   }

   // Now process the context TODO 086
   double *dAngleOut = new double[m_ulNumCells];
   nRet = GDALGridContextProcess(pContext, 0, m_nXGridSize - 1, 0, m_nYGridSize - 1, m_nXGridSize, m_nYGridSize, GDT_Float64, dAngleOut, NULL, NULL);

   if (nRet == CE_Failure)
   {
      delete[] dHeightOut;
      delete[] dAngleOut;
      delete pOptions;
      return RTN_ERR_GRIDCREATE;
   }

   // Get rid of the context
   GDALGridContextFree(pContext);

   // OK, now create a gridded version of wave period: first create the GDAL context
   // pContext = GDALGridContextCreate(GGA_InverseDistanceToAPower, pOptions, nUserPoints, &m_VdDeepWaterWaveStationX[0], &m_VdDeepWaterWaveStationY[0], &m_VdThisIterDeepWaterWaveStationPeriod[0], true);
   pContext = GDALGridContextCreate(GGA_InverseDistanceToAPower, pOptions, nUserPoints, m_VdDeepWaterWaveStationX.data(), m_VdDeepWaterWaveStationY.data(), m_VdThisIterDeepWaterWaveStationPeriod.data(), true);

   if (pContext == NULL)
   {
      delete pOptions;
      return RTN_ERR_GRIDCREATE;
   }

   // Now process the context TODO 086
   double *dPeriopdOut = new double[m_ulNumCells];
   nRet = GDALGridContextProcess(pContext, 0, m_nXGridSize - 1, 0, m_nYGridSize - 1, m_nXGridSize, m_nYGridSize, GDT_Float64, dPeriopdOut, NULL, NULL);

   if (nRet == CE_Failure)
   {
      delete[] dPeriopdOut;
      delete pOptions;
      return RTN_ERR_GRIDCREATE;
   }

   // Get rid of the context
   GDALGridContextFree(pContext);

   // The output from GDALGridCreate() is in dHeightOut, dAngleOut and dPeriopdOut but must be reversed
   vector<double> VdHeight;
   vector<double> VdAngle;
   vector<double> VdPeriod;

   int n = 0;
   int nValidHeight = 0;
   int nValidAngle = 0;
   int nValidPeriod = 0;

   double dAvgHeight = 0;
   double dAvgAngle = 0;
   double dAvgPeriod = 0;

   for (int nY = m_nYGridSize - 1; nY >= 0; nY--)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         if (isfinite(dHeightOut[n]))
         {
            VdHeight.push_back(dHeightOut[n]);

            dAvgHeight += dHeightOut[n];
            nValidHeight++;
         }

         else
         {
            VdHeight.push_back(m_dMissingValue);
         }

         if (isfinite(dAngleOut[n]))
         {
            VdAngle.push_back(dAngleOut[n]);

            dAvgAngle += dAngleOut[n];
            nValidAngle++;
         }

         else
         {
            VdAngle.push_back(m_dMissingValue);
         }

         if (isfinite(dPeriopdOut[n]))
         {
            VdPeriod.push_back(dPeriopdOut[n]);

            dAvgPeriod += dPeriopdOut[n];
            nValidPeriod++;
         }

         else
         {
            VdPeriod.push_back(m_dMissingValue);
         }

         // LogStream << " nX = " << nX << " nY = " << nY << " n = " << n << " dHeightOut[n] = " << dHeightOut[n] << " dAngleOut[n] = " << dAngleOut[n] << endl;
         n++;
      }
   }

   // Calculate averages
   dAvgHeight /= nValidHeight;
   dAvgAngle /= nValidAngle;
   dAvgPeriod /= nValidPeriod;

   // Tidy
   delete pOptions;
   delete[] dHeightOut;
   delete[] dAngleOut;
   delete[] dPeriopdOut;

   // Now update all raster cells
   n = 0;

   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         if (bFPIsEqual(VdHeight[n], m_dMissingValue, TOLERANCE))
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWaveHeight(dAvgHeight);
         else
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWaveHeight(VdHeight[n]);

         if (bFPIsEqual(VdAngle[n], m_dMissingValue, TOLERANCE))
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWaveAngle(dAvgAngle);
         else
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWaveAngle(
               VdAngle[n]);

         if (bFPIsEqual(VdPeriod[n], m_dMissingValue, TOLERANCE))
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWavePeriod(dAvgPeriod);
         else
            m_pRasterGrid->m_Cell[nX][nY].SetCellDeepWaterWavePeriod(
               VdPeriod[n]);

         // LogStream << " [" << nX << "][" << nY << "] deep water wave height = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight() << " deep water wave angle = " << m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle() << endl;
         n++;
      }
   }

   // // DEBUG CODE
   // ===========================================================================================================
   // string strOutFile = m_strOutPath;
   // strOutFile += "init_deep_water_wave_height_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   // GDALDriver* pDriver =
   // GetGDALDriverManager()->GetDriverByName("gtiff"); GDALDataset*
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize,
   // m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_ulNumCells];
   // int nn = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   //          // Write this value to the array
   // pdRaster[nn] =
   // m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight(); nn++;
   // }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_nMissingValue);
   // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize,
   // pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   // // DEBUG CODE
   // ===========================================================================================================

   // // DEBUG CODE
   // ===========================================================================================================
   // strOutFile = m_strOutPath;
   // strOutFile += "init_deep_water_wave_angle_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize,
   // m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // nn = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   //          // Write this value to the array
   // pdRaster[nn] =
   // m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle(); nn++;
   // }
   // }
   //
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_nMissingValue);
   // nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize,
   // pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // return RTN_ERR_GRIDCREATE;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE
   // ===========================================================================================================

   return RTN_OK;
}
