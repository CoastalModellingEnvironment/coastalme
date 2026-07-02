/*!
   \file locate_coast.cpp
   \brief Finds the coastline on the raster grid
   \details TODO 001 A more detailed description of these routines.
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
#include <assert.h>

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

// #include <ios>
// using std::fixed;

// #include <string>
// using std::to_string;

// #include <iomanip>
// using std::setprecision;

#include <stack>
using std::stack;

// #include <cpl_conv.h>
// #include <cpl_error.h>
// #include <cpl_string.h>
// #include <gdal.h>
// #include <gdal_alg.h>
// #include <gdal_priv.h>

#include "cme.h"
#include "2di_point.h"
#include "i_line.h"
#include "line.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! First find all connected sea areas, then locate the vector coastline(s), then put these onto the raster grid
//===============================================================================================================================
int CSimulation::nLocateSeaAndCoasts(void)
{
   // Find all connected sea cells
   FindAllSeaCells();

   // Find every coastline on the raster grid, mark raster cells, then create the vector coastline
   int const nRet = nTraceAllCoasts();
   if (nRet != RTN_OK)
      return nRet;

   // Have we created any coasts?
   if (m_VCoast.empty())
   {
      cerr << m_ulIter << ": " << ERR << "no coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;

      return RTN_ERR_NO_COAST;
   }

   // Is this the highest SWL so far? If so, save this for all coasts
   if (m_bHighestSWLSoFar)
   {
      m_VHighestSWLCoastLine.clear();

      for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      {
         CGeomLine LCoast;

         LCoast = *m_VCoast[nCoast].pLGetCoastlineExtCRS();
         m_VHighestSWLCoastLine.push_back(LCoast);
      }
   }

   // Is this the lowest SWL so far? If so, save this for all coasts
   if (m_bLowestSWLSoFar)
   {
      m_VLowestSWLCoastLine.clear();

      for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      {
         CGeomLine LCoast;

         LCoast = *m_VCoast[nCoast].pLGetCoastlineExtCRS();
         m_VLowestSWLCoastLine.push_back(LCoast);
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Finds and flags all sea areas which have at least one cell at a grid edge (i.e. does not flag 'inland' seas)
//===============================================================================================================================
void CSimulation::FindAllSeaCells(void)
{
   // Go along each list of edge cells, north edge first
   if (! m_bOmitSearchNorthEdge)
   {
      if (m_bSearchNorthEdgeForward)
      {
         for (int n = 0; n < static_cast<int>(m_VPtiNorthEdgeCell.size()); n++)
         {
            int const nX = m_VPtiNorthEdgeCell[n].nGetX();
            int const nY = m_VPtiNorthEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
      else
      {
         for (int n = static_cast<int>(m_VPtiNorthEdgeCell.size())-1; n >= 0; n--)
         {
            int const nX = m_VPtiNorthEdgeCell[n].nGetX();
            int const nY = m_VPtiNorthEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
   }

   // Now go along the south edge cells
   if (! m_bOmitSearchSouthEdge)
   {
      if (m_bSearchSouthEdgeForward)
      {
         for (int n = 0; n < static_cast<int>(m_VPtiSouthEdgeCell.size()); n++)
         {
            int const nX = m_VPtiSouthEdgeCell[n].nGetX();
            int const nY = m_VPtiSouthEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
      else
      {
         for (int n = static_cast<int>(m_VPtiSouthEdgeCell.size())-1; n >= 0; n--)
         {
            int const nX = m_VPtiSouthEdgeCell[n].nGetX();
            int const nY = m_VPtiSouthEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
   }

   // Now go along the west edge cells
   if (! m_bOmitSearchWestEdge)
   {
      if (m_bSearchWestEdgeForward)
      {
         for (int n = 0; n < static_cast<int>(m_VPtiWestEdgeCell.size()); n++)
         {
            int const nX = m_VPtiWestEdgeCell[n].nGetX();
            int const nY = m_VPtiWestEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
      else
      {
         for (int n = static_cast<int>(m_VPtiWestEdgeCell.size())-1; n >= 0; n--)
         {
            int const nX = m_VPtiWestEdgeCell[n].nGetX();
            int const nY = m_VPtiWestEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
   }

   // Finally go along the east edge cells
   if (! m_bOmitSearchEastEdge)
   {
      if (m_bSearchEastEdgeForward)
      {
         for (int n = 0; n < static_cast<int>(m_VPtiEastEdgeCell.size()); n++)
         {
            int const nX = m_VPtiEastEdgeCell[n].nGetX();
            int const nY = m_VPtiEastEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
      else
      {
         for (int n = static_cast<int>(m_VPtiEastEdgeCell.size())-1; n >= 0; n--)
         {
            int const nX = m_VPtiEastEdgeCell[n].nGetX();
            int const nY = m_VPtiEastEdgeCell[n].nGetY();

            if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
            {
               // This edge cell is below SWL and sea depth remains set to zero
               CellByCellFillSea(nX, nY);

               return;
            }
         }
      }
   }
}

//===============================================================================================================================
//! Cell-by-cell fills all sea cells starting from a given cell. The cell-by-cell fill (aka 'floodfill') code used here is adapted from an example by Lode Vandevenne (http://lodev.org/cgtutor/floodfill.html#Scanline_Floodfill_Algorithm_With_Stack)
//===============================================================================================================================
void CSimulation::CellByCellFillSea(int const nXStart, int const nYStart)
{
   // For safety check
   int const nRoundLoopMax = m_nXGridSize * m_nYGridSize;

   // Create an empty stack
   stack<CGeom2DIPoint> PtiStack;

   // Start at the given edge cell, push this onto the stack
   PtiStack.push(CGeom2DIPoint(nXStart, nYStart));

   // Then do the cell-by-cell fill loop until there are no more cell coordinates on the stack
   int nRoundLoop = 0;

   while (! PtiStack.empty())
   {
      // Safety check
      if (nRoundLoop++ > nRoundLoopMax)
         break;

      CGeom2DIPoint const Pti = PtiStack.top();
      PtiStack.pop();

      int nX = Pti.nGetX();
      int const nY = Pti.nGetY();

      while ((nX >= 0) && (!m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()))
         nX--;

      nX++;

      bool bSpanAbove = false;
      bool bSpanBelow = false;

      while ((nX < m_nXGridSize) && (!m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
      {
         // Set the sea depth for this cell
         m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

         CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetCellLandform();
         int const nCat = pLandform->nGetLandformCategory();

         // Have we had sediment input here?
         if ((nCat == LF_SEDIMENT_INPUT_CONSOLIDATED) || (nCat == LF_SEDIMENT_INPUT_UNCONSOLIDATED))
         {
            if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
            {
               m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();

               // Set this sea cell to have deep water (off-shore) wave orientation and height, will change this later for cells closer to the shoreline if we have on-shore waves
               m_pRasterGrid->m_Cell[nX][nY].SetWaveValuesToDeepWaterWaveValues();
            }
         }
         else
         {
            // No sediment input here, just mark as sea
            m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();
            pLandform->SetLandformCategory(LF_SEA);

            // Set this sea cell to have deep water (off-shore) wave orientation and height, will change this later for cells closer to the shoreline if we have on-shore waves
            m_pRasterGrid->m_Cell[nX][nY].SetWaveValuesToDeepWaterWaveValues();
         }

         // Now sort out the x-y extremities of the contiguous sea for the bounding box (used later in wave propagation)
         if (nX < m_nXMinBoundingBox)
            m_nXMinBoundingBox = nX;

         if (nX > m_nXMaxBoundingBox)
            m_nXMaxBoundingBox = nX;

         if (nY < m_nYMinBoundingBox)
            m_nYMinBoundingBox = nY;

         if (nY > m_nYMaxBoundingBox)
            m_nYMaxBoundingBox = nY;

         // Update count
         m_ulThisIterNumSeaCells++;

         if ((! bSpanAbove) && (nY > 0) && (!m_pRasterGrid->m_Cell[nX][nY - 1].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY - 1].bIsInundated()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY - 1));
            bSpanAbove = true;
         }
         else if (bSpanAbove && (nY > 0) && (!m_pRasterGrid->m_Cell[nX][nY - 1].bBasementElevIsMissingValue()) && (!m_pRasterGrid->m_Cell[nX][nY - 1].bIsInundated()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && (nY < m_nYGridSize - 1) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY + 1].bIsInundated()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY + 1));
            bSpanBelow = true;
         }
         else if (bSpanBelow && (nY < m_nYGridSize - 1) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bBasementElevIsMissingValue()) && (!m_pRasterGrid->m_Cell[nX][nY + 1].bIsInundated()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath + "is_contiguous_sea_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   //
   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   //    for (int nX = 0; nX < m_nXGridSize; nX++)
   //    {
   //    pdRaster[n++] = m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea();
   //    }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   // return;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath + "is_inundated_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   //
   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
   //
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   // pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   //    for (int nX = 0; nX < m_nXGridSize; nX++)
   //    {
   //       pdRaster[n++] = m_pRasterGrid->m_Cell[nX][nY].bIsInundated();
   //    }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   // return;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // LogStream << m_ulIter << ": cell-by-cell fill of sea from [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "} with SWL = " << m_dThisIterSWL << ", " << m_ulThisIterNumSeaCells << " of " << m_ulNumCells << " cells now marked as sea (" <<  fixed << setprecision(3) << 100.0 * m_ulThisIterNumSeaCells / m_ulNumCells << " %)" << endl;
   //
   // LogStream << " m_nXMinBoundingBox = " << m_nXMinBoundingBox << " m_nXMaxBoundingBox = " << m_nXMaxBoundingBox << " m_nYMinBoundingBox = " << m_nYMinBoundingBox << " m_nYMaxBoundingBox = " << m_nYMaxBoundingBox << endl;
   // // DEBUG CODE ===========================================================================================================
}

//===============================================================================================================================
//! Locates all the potential coastline start/finish points on the edges of the raster grid, then traces vector coastline(s) from these start points
//===============================================================================================================================
int CSimulation::nTraceAllCoasts(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Tracing coasts" << endl;

   int nValidCoast = 0;

   vector<CGeom2DIPoint> V2DIPossibleStartCell;
   vector<int> VnPossibleStartCellEdge;
   vector<int> VnPossibleStartCellHandedness;
   vector<int> VnPossibleStartCellSearchDirection;
   vector<bool> VbTraced;

   // Go along each list of edge cells from low to high elevation, so that the most seaward possible coast point is found first. Start with north edge
   if (! m_bOmitSearchNorthEdge)
   {
      if (m_bSearchNorthEdgeForward)
      {
         // Searching forward (i.e. in direction of increasing indices)
         for (int n = 0; n < static_cast<int>(m_VPtiNorthEdgeCell.size())-1; n++)
         {
            int const nXThis = m_VPtiNorthEdgeCell[n].nGetX();
            int const nYThis = m_VPtiNorthEdgeCell[n].nGetY();
            int const nXNext = m_VPtiNorthEdgeCell[n + 1].nGetX();
            int const nYNext = m_VPtiNorthEdgeCell[n + 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(NORTH);
               VnPossibleStartCellHandedness.push_back(RIGHT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(SOUTH);
               VbTraced.push_back(false);
            }
         }
      }
      else
      {
         // Searching backward (i.e. in direction of decreasing indices)
         for (int n = static_cast<int>(m_VPtiNorthEdgeCell.size())-1; n > 0; n--)
         {
            int const nXThis = m_VPtiNorthEdgeCell[n].nGetX();
            int const nYThis = m_VPtiNorthEdgeCell[n].nGetY();
            int const nXNext = m_VPtiNorthEdgeCell[n - 1].nGetX();
            int const nYNext = m_VPtiNorthEdgeCell[n - 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(NORTH);
               VnPossibleStartCellHandedness.push_back(LEFT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(SOUTH);
               VbTraced.push_back(false);
            }
         }
      }
   }

   // Now go along south edge
   if (! m_bOmitSearchSouthEdge)
   {
      if (m_bSearchSouthEdgeForward)
      {
         // Searching forward (i.e. in direction of increasing indices)
         for (int n = 0; n < static_cast<int>(m_VPtiSouthEdgeCell.size())-1; n++)
         {
            int const nXThis = m_VPtiSouthEdgeCell[n].nGetX();
            int const nYThis = m_VPtiSouthEdgeCell[n].nGetY();
            int const nXNext = m_VPtiSouthEdgeCell[n + 1].nGetX();
            int const nYNext = m_VPtiSouthEdgeCell[n + 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(SOUTH);
               VnPossibleStartCellHandedness.push_back(LEFT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(NORTH);
               VbTraced.push_back(false);
            }
         }
      }
      else
      {
         // Searching backward (i.e. in direction of decreasing indices)
         for (int n = static_cast<int>(m_VPtiSouthEdgeCell.size())-1; n > 0; n--)
         {
            int const nXThis = m_VPtiSouthEdgeCell[n].nGetX();
            int const nYThis = m_VPtiSouthEdgeCell[n].nGetY();
            int const nXNext = m_VPtiSouthEdgeCell[n - 1].nGetX();
            int const nYNext = m_VPtiSouthEdgeCell[n - 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(SOUTH);
               VnPossibleStartCellHandedness.push_back(RIGHT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(NORTH);
               VbTraced.push_back(false);
            }
         }
      }
   }

   // Now go along west edge
   if (! m_bOmitSearchWestEdge)
   {
      if (m_bSearchWestEdgeForward)
      {
         // Searching forward (i.e. in direction of increasing indices)
         for (int n = 0; n < static_cast<int>(m_VPtiWestEdgeCell.size())-1; n++)
         {
            int const nXThis = m_VPtiWestEdgeCell[n].nGetX();
            int const nYThis = m_VPtiWestEdgeCell[n].nGetY();
            int const nXNext = m_VPtiWestEdgeCell[n + 1].nGetX();
            int const nYNext = m_VPtiWestEdgeCell[n + 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(WEST);
               VnPossibleStartCellHandedness.push_back(LEFT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(EAST);
               VbTraced.push_back(false);
            }
         }
      }
      else
      {
         // Searching backward (i.e. in direction of decreasing indices)
         for (int n = static_cast<int>(m_VPtiWestEdgeCell.size())-1; n > 0; n--)
         {
            int const nXThis = m_VPtiWestEdgeCell[n].nGetX();
            int const nYThis = m_VPtiWestEdgeCell[n].nGetY();
            int const nXNext = m_VPtiWestEdgeCell[n - 1].nGetX();
            int const nYNext = m_VPtiWestEdgeCell[n - 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(WEST);
               VnPossibleStartCellHandedness.push_back(RIGHT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(EAST);
               VbTraced.push_back(false);
            }
         }
      }
   }

   // Finally go along east edge
   if (! m_bOmitSearchEastEdge)
   {
      if (m_bSearchEastEdgeForward)
      {
         // Searching forward (i.e. in direction of increasing indices)
         for (int n = 0; n < static_cast<int>(m_VPtiEastEdgeCell.size())-1; n++)
         {
            int const nXThis = m_VPtiEastEdgeCell[n].nGetX();
            int const nYThis = m_VPtiEastEdgeCell[n].nGetY();
            int const nXNext = m_VPtiEastEdgeCell[n + 1].nGetX();
            int const nYNext = m_VPtiEastEdgeCell[n + 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(EAST);
               VnPossibleStartCellHandedness.push_back(RIGHT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(WEST);
               VbTraced.push_back(false);
            }
         }
      }
      else
      {
         // Searching backward (i.e. in direction of decreasing indices)
         for (int n = static_cast<int>(m_VPtiEastEdgeCell.size())-1; n > 0; n--)
         {
            int const nXThis = m_VPtiEastEdgeCell[n].nGetX();
            int const nYThis = m_VPtiEastEdgeCell[n].nGetY();
            int const nXNext = m_VPtiEastEdgeCell[n - 1].nGetX();
            int const nYNext = m_VPtiEastEdgeCell[n - 1].nGetY();

            if (bIdentifyPossibleCoastStart(nXThis, nYThis, nXNext, nYNext, &V2DIPossibleStartCell))
            {
               VnPossibleStartCellEdge.push_back(EAST);
               VnPossibleStartCellHandedness.push_back(LEFT_HANDED);
               VnPossibleStartCellSearchDirection.push_back(WEST);
               VbTraced.push_back(false);
            }
         }
      }
   }

   // Any possible coastline start/finish cells found?
   if (V2DIPossibleStartCell.size() == 0)
   {
      LogStream << m_ulIter << ": no coastline start/finish points found after grid edges searched.";

      if (m_bOmitSearchNorthEdge || m_bOmitSearchSouthEdge || m_bOmitSearchWestEdge || m_bOmitSearchEastEdge)
      {
         LogStream << " Note that the following grid edges were not searched: " << (m_bOmitSearchNorthEdge ? "N " : "") << (m_bOmitSearchSouthEdge ? "S " : "") << (m_bOmitSearchWestEdge ? "W " : "") << (m_bOmitSearchEastEdge ? "E " : "");
      }

      LogStream << endl;

      return RTN_ERR_NO_START_FINISH_POINTS_TRACING_COAST;
   }


   // All OK, now trace from each of these possible start/finish points
   for (int n = 0; n < static_cast<int>(V2DIPossibleStartCell.size()); n++)
   {
      if (! VbTraced[n])
      {

         int const nRet = nTraceCoastLine(n, &V2DIPossibleStartCell, &VbTraced, VnPossibleStartCellEdge[n], VnPossibleStartCellHandedness[n], VnPossibleStartCellSearchDirection[n]);
         if (nRet == RTN_OK)
         {
            // We have a valid coastline starting from this possible start cell
            nValidCoast++;
         }
      }
   }

   if (nValidCoast == 0)
   {
      // Still no valid coasts found, so we have to give up
      cerr << m_ulIter << ": no valid coasts found, see " << m_strLogFile << " for more information" << endl;
      return RTN_ERR_NO_VALID_COAST;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Identifies a possible start- or end-of-coast edge cell, searching from low elevation (sea) towards high elevation (land)
//===============================================================================================================================
bool CSimulation::bIdentifyPossibleCoastStart(int const nXThis, int const nYThis, int const nXNext, int const nYNext, vector<CGeom2DIPoint>* pV2DIPossibleStartCell)
{
   // Get "Is it sea?" information for 'this' and 'next' cells
   bool const bThisCellIsSea = m_pRasterGrid->m_Cell[nXThis][nYThis].bIsInContiguousSea();
   bool const bNextCellIsSea = m_pRasterGrid->m_Cell[nXNext][nYNext].bIsInContiguousSea();

   // We are searching from sea to land: so is this cell sea and the next cell land?
   if (bThisCellIsSea && (! bNextCellIsSea))
   {
      // All OK, so flag the 'next' cell
      m_pRasterGrid->m_Cell[nXNext][nYNext].SetPossibleCoastStartCell();

      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t flagging [" << nXNext << "][" << nYNext << "] = {" << dGridCentroidXToExtCRSX(nXNext) << ", " << dGridCentroidYToExtCRSY(nYNext) << "} as possible coast start or end cell" << endl;

      // And save it as a possible coastline start/end cell
      pV2DIPossibleStartCell->push_back(CGeom2DIPoint(nXNext, nYNext));

      return true;
   }

   return false;
}

//===============================================================================================================================
//! Traces a coastline (which is defined to be just above still water level) on the grid using the 'wall follower' rule for maze traversal (http://en.wikipedia.org/wiki/Maze_solving_algorithm#Wall_follower). The resulting vector coastline is then smoothed
//===============================================================================================================================
int CSimulation::nTraceCoastLine(int const nTraceFromStartCellIndex, vector<CGeom2DIPoint> const* pV2DIPossibleStartCell, vector<bool>* pVbTraced, int const nStartEdge, int const nHandedness, int const nStartSearchDirection)
{
   // bool bHitStartCell = false;
   // bool bOnCoast = false;
   bool bHasLeftStartEdge = false;
   // bool bTooLong = false;
   bool bDeadEnd = false;
   bool const bHereBefore = false;

   int const nStartX = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetX();
   int const nStartY = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetY();
   int nX = nStartX;
   int nY = nStartY;
   int nSearchDirection = nStartSearchDirection;
   int nRoundLoop = -1;

   // Temporary coastline as integer points (grid CRS)
   CGeomILine ILTempGridCRS;

   // Start at this grid-edge point and trace the rest of the coastline using the 'wall follower' rule for maze traversal, trying to keep next to cells flagged as sea
   do
   {
      nRoundLoop++;

      // Have we hit a grid edge?
      if (m_pRasterGrid->m_Cell[nX][nY].bIsBoundingBoxEdge())
      {
         // We have, but is a different edge from the start edge?
         int const nEdge = m_pRasterGrid->m_Cell[nX][nY].nGetBoundingBoxEdge();
         if (nEdge != nStartEdge)
         {
            // We have hit a different edge, so this is a possible coast
            break;
         }
      }

      // Safety check: have we returned to the start point?
      if ((nRoundLoop > 1) && (nX == nStartX) && (nY == nStartY))
      {
         // We have returned to the start point
         break;
      }

      // // Another safety check: have we gone on too long?
      // if (nRoundLoop > m_nCoastMax)
      // {
      //    bTooLong = true;
      //
      //    // LogStream << m_ulIter << ":\t abandoning possible coastline, traced from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, exceeded maximum search length (" << m_nCoastMax << ")" << endl;
      //
      //    // for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
      //    // LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
      //    // LogStream << endl;
      //
      //    break;
      // }

      // Another safety check: have we visited this cell before?
      // if (ILTempGridCRS.bIsPresent(nX, nY))
      // {
      //    // We've been here before
      //    bHereBefore = true;
      //
      //    LogStream << m_ulIter << ":\t abandoning possible coastline, traced from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, since have already visited [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
      //
      //    break;
      // }

      // OK so far: so have we left the start edge?
      if (! bHasLeftStartEdge)
      {
         if (((nStartSearchDirection == SOUTH) && (nY > nStartY)) || ((nStartSearchDirection == NORTH) && (nY < nStartY)) || ((nStartSearchDirection == EAST) && (nX > nStartX)) || ((nStartSearchDirection == WEST) && (nX < nStartX)))
         {
            // LogStream << "Left start edge with nX = " << nX << " nY = " << nY << endl;
            bHasLeftStartEdge = true;
         }
      }

      // LogStream << "Now at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} bHasLeftStartEdge = " << bHasLeftStartEdge << " bOnCoast = " << bOnCoast << endl;

      // // If the vector coastline has left the start edge, and we hit a possible coast start point from which a coastline has not yet been traced, then hooray! We've traced a coast. Leave the loop
      // if (bHasLeftStartEdge && bOnCoast)
      // {
      //    for (int nn = 0; nn < static_cast<int>(pVbTraced->size()); nn++)
      //    {
      //       bool const bTraced = pVbTraced->at(nn);
      //       if ((nn != nTraceFromStartCellIndex) && (! bTraced))
      //       {
      //          int const nXPoss = pV2DIPossibleStartCell->at(nn).nGetX();
      //          int const nYPoss = pV2DIPossibleStartCell->at(nn).nGetY();
      //
      //          LogStream << "In 'Leave the edge' loop for [" << nX << "][" << nY << "] bTraced = " << bTraced << " nn = " << nn << " nTraceFromStartCellIndex = " << nTraceFromStartCellIndex << " possible start cell = [" << nXPoss << "][" << nYPoss << "]" << endl;
      //
      //          if ((nX == nXPoss) && (nY == nYPoss))
      //          {
      //             if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      //                LogStream << m_ulIter << ":\t possible coastline found, traced from [" << nStartX << "][" << nStartY << "]  = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, ended at possible coast start cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
      //
      //             pVbTraced->at(nn) = true;
      //             bHitStartCell = true;
      //             break;
      //          }
      //       }
      //    }
      // }
      //
      // if (bHitStartCell)
      //    break;

      // OK now prepare for the next iteration of the search
      int nXSeaward = 0;
      int nYSeaward = 0;
      int nSeawardNewDirection = 0;
      int nXStraightOn = 0;
      int nYStraightOn = 0;
      int nXAntiSeaward = 0;
      int nYAntiSeaward = 0;
      int nAntiSeawardNewDirection = 0;
      int nXGoBack = 0;
      int nYGoBack = 0;
      int nGoBackNewDirection = 0;

      CGeom2DIPoint const Pti(nX, nY);

      // Set up the variables
      switch (nHandedness)
      {
      case RIGHT_HANDED:
         // The sea is to the right-hand side of the coast as we traverse it. We are just inland, so we need to keep heading right to find the sea
         switch (nSearchDirection)
         {
         case NORTH:
            // The sea is towards the RHS (E) of the coast, so first try to go right (to the E)
            nXSeaward = nX + 1;
            nYSeaward = nY;
            nSeawardNewDirection = EAST;

            // If can't do this, try to go straight on (to the N)
            nXStraightOn = nX;
            nYStraightOn = nY - 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (W)
            nXAntiSeaward = nX - 1;
            nYAntiSeaward = nY;
            nAntiSeawardNewDirection = WEST;

            // As a last resort, go back (to the S)
            nXGoBack = nX;
            nYGoBack = nY + 1;
            nGoBackNewDirection = SOUTH;

            break;

         case NORTH_EAST:
            // The sea is towards the RHS (SE) of the coast, so first try to go right (to the SE)
            nXSeaward = nX + 1;
            nYSeaward = nY + 1;
            nSeawardNewDirection = SOUTH_EAST;

            // If can't do this, try to go straight on (to the NE)
            nXStraightOn = nX + 1;
            nYStraightOn = nY - 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (NW)
            nXAntiSeaward = nX - 1;
            nYAntiSeaward = nY - 1;
            nAntiSeawardNewDirection = NORTH_WEST;

            // As a last resort, go back (to the SW)
            nXGoBack = nX - 1;
            nYGoBack = nY + 1;
            nGoBackNewDirection = SOUTH_WEST;

            break;

         case EAST:
            // The sea is towards the RHS (S) of the coast, so first try to go right (to the S)
            nXSeaward = nX;
            nYSeaward = nY + 1;
            nSeawardNewDirection = SOUTH;

            // If can't do this, try to go straight on (to the E)
            nXStraightOn = nX + 1;
            nYStraightOn = nY;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (N)
            nXAntiSeaward = nX;
            nYAntiSeaward = nY - 1;
            nAntiSeawardNewDirection = NORTH;

            // As a last resort, go back (to the W)
            nXGoBack = nX - 1;
            nYGoBack = nY;
            nGoBackNewDirection = WEST;

            break;

         case SOUTH_EAST:
            // The sea is towards the RHS (SW) of the coast, so first try to go right (to the SW)
            nXSeaward = nX - 1;
            nYSeaward = nY + 1;
            nSeawardNewDirection = SOUTH_WEST;

            // If can't do this, try to go straight on (to the SE)
            nXStraightOn = nX + 1;
            nYStraightOn = nY + 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (NE)
            nXAntiSeaward = nX + 1;
            nYAntiSeaward = nY - 1;
            nAntiSeawardNewDirection = NORTH_EAST;

            // As a last resort, go back (to the NW)
            nXGoBack = nX - 1;
            nYGoBack = nY - 1;
            nGoBackNewDirection = NORTH_WEST;

            break;

         case SOUTH:
            // The sea is towards the RHS (W) of the coast, so first try to go right (to the W)
            nXSeaward = nX - 1;
            nYSeaward = nY;
            nSeawardNewDirection = WEST;

            // If can't do this, try to go straight on (to the S)
            nXStraightOn = nX;
            nYStraightOn = nY + 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (E)
            nXAntiSeaward = nX + 1;
            nYAntiSeaward = nY;
            nAntiSeawardNewDirection = EAST;

            // As a last resort, go back (to the N)
            nXGoBack = nX;
            nYGoBack = nY - 1;
            nGoBackNewDirection = NORTH;

            break;

         case SOUTH_WEST:
            // The sea is towards the RHS (NW) of the coast, so first try to go right (to the NW)
            nXSeaward = nX - 1;
            nYSeaward = nY - 1;
            nSeawardNewDirection = NORTH_WEST;

            // If can't do this, try to go straight on (to the SW)
            nXStraightOn = nX - 1;
            nYStraightOn = nY + 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (SE)
            nXAntiSeaward = nX + 1;
            nYAntiSeaward = nY + 1;
            nAntiSeawardNewDirection = SOUTH_EAST;

            // As a last resort, go back (to the NE)
            nXGoBack = nX + 1;
            nYGoBack = nY - 1;
            nGoBackNewDirection = NORTH_EAST;

            break;

         case WEST:
            // The sea is towards the RHS (N) of the coast, so first try to go right (to the N)
            nXSeaward = nX;
            nYSeaward = nY - 1;
            nSeawardNewDirection = NORTH;

            // If can't do this, try to go straight on (to the W)
            nXStraightOn = nX - 1;
            nYStraightOn = nY;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (S)
            nXAntiSeaward = nX;
            nYAntiSeaward = nY + 1;
            nAntiSeawardNewDirection = SOUTH;

            // As a last resort, go back (to the E)
            nXGoBack = nX + 1;
            nYGoBack = nY;
            nGoBackNewDirection = EAST;

            break;

         case NORTH_WEST:
            // The sea is towards the RHS (NE) of the coast, so first try to go right (to the NE)
            nXSeaward = nX + 1;
            nYSeaward = nY - 1;
            nSeawardNewDirection = NORTH_EAST;

            // If can't do this, try to go straight on (to the NW)
            nXStraightOn = nX - 1;
            nYStraightOn = nY - 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the LHS (SW)
            nXAntiSeaward = nX - 1;
            nYAntiSeaward = nY + 1;
            nAntiSeawardNewDirection = SOUTH_WEST;

            // As a last resort, go back (to the SE)
            nXGoBack = nX + 1;
            nYGoBack = nY - 1;
            nGoBackNewDirection = SOUTH_EAST;

            break;
         }

         break;

      case LEFT_HANDED:

         // The sea is to the left-hand side of the coast as we traverse it. We are just inland, so we need to keep heading left to find the sea
         switch (nSearchDirection)
         {
         case NORTH:
            // The sea is towards the LHS (W) of the coast, so first try to go left (to the W)
            nXSeaward = nX - 1;
            nYSeaward = nY;
            nSeawardNewDirection = WEST;

            // If can't do this, try to go straight on (to the N)
            nXStraightOn = nX;
            nYStraightOn = nY - 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (E)
            nXAntiSeaward = nX + 1;
            nYAntiSeaward = nY;
            nAntiSeawardNewDirection = EAST;

            // As a last resort, go back (to the S)
            nXGoBack = nX;
            nYGoBack = nY + 1;
            nGoBackNewDirection = SOUTH;

            break;

         case NORTH_EAST:
            // The sea is towards the LHS (NW) of the coast, so first try to go left (to the NW)
            nXSeaward = nX - 1;
            nYSeaward = nY - 1;
            nSeawardNewDirection = NORTH_WEST;

            // If can't do this, try to go straight on (to the NE)
            nXStraightOn = nX + 1;
            nYStraightOn = nY - 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (SE)
            nXAntiSeaward = nX + 1;
            nYAntiSeaward = nY + 1;
            nAntiSeawardNewDirection = SOUTH_EAST;

            // As a last resort, go back (to the SW)
            nXGoBack = nX - 1;
            nYGoBack = nY + 1;
            nGoBackNewDirection = SOUTH_WEST;

            break;

         case EAST:
            // The sea is towards the LHS (N) of the coast, so first try to go left (to the N)
            nXSeaward = nX;
            nYSeaward = nY - 1;
            nSeawardNewDirection = NORTH;

            // If can't do this, try to go straight on (to the E)
            nXStraightOn = nX + 1;
            nYStraightOn = nY;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (S)
            nXAntiSeaward = nX;
            nYAntiSeaward = nY + 1;
            nAntiSeawardNewDirection = SOUTH;

            // As a last resort, go back (to the W)
            nXGoBack = nX - 1;
            nYGoBack = nY;
            nGoBackNewDirection = WEST;

            break;

         case SOUTH_EAST:
            // The sea is towards the LHS (NE) of the coast, so first try to go left (to the NE)
            nXSeaward = nX + 1;
            nYSeaward = nY - 1;
            nSeawardNewDirection = NORTH_EAST;

            // If can't do this, try to go straight on (to the SE)
            nXStraightOn = nX + 1;
            nYStraightOn = nY + 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (SW)
            nXAntiSeaward = nX - 1;
            nYAntiSeaward = nY + 1;
            nAntiSeawardNewDirection = SOUTH_WEST;

            // As a last resort, go back (to the NW)
            nXGoBack = nX - 1;
            nYGoBack = nY - 1;
            nGoBackNewDirection = NORTH_WEST;

            break;

         case SOUTH:
            // The sea is towards the LHS (E) of the coast, so first try to go left (to the E)
            nXSeaward = nX + 1;
            nYSeaward = nY;
            nSeawardNewDirection = EAST;

            // If can't do this, try to go straight on (to the S)
            nXStraightOn = nX;
            nYStraightOn = nY + 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (W)
            nXAntiSeaward = nX - 1;
            nYAntiSeaward = nY;
            nAntiSeawardNewDirection = WEST;

            // As a last resort, go back (to the N)
            nXGoBack = nX;
            nYGoBack = nY - 1;
            nGoBackNewDirection = NORTH;

            break;

         case SOUTH_WEST:
            // The sea is towards the LHS (SE) of the coast, so first try to go left (to the SE)
            nXSeaward = nX + 1;
            nYSeaward = nY + 1;
            nSeawardNewDirection = SOUTH_EAST;

            // If can't do this, try to go straight on (to the SW)
            nXStraightOn = nX - 1;
            nYStraightOn = nY + 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (NW)
            nXAntiSeaward = nX - 1;
            nYAntiSeaward = nY - 1;
            nAntiSeawardNewDirection = NORTH_WEST;

            // As a last resort, go back (to the NE)
            nXGoBack = nX + 1;
            nYGoBack = nY - 1;
            nGoBackNewDirection = NORTH_EAST;

            break;

         case WEST:
            // The sea is towards the LHS (S) of the coast, so first try to go left (to the S)
            nXSeaward = nX;
            nYSeaward = nY + 1;
            nSeawardNewDirection = SOUTH;

            // If can't do this, try to go straight on (to the W)
            nXStraightOn = nX - 1;
            nYStraightOn = nY;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (N)
            nXAntiSeaward = nX;
            nYAntiSeaward = nY - 1;
            nAntiSeawardNewDirection = NORTH;

            // As a last resort, go back (to the E)
            nXGoBack = nX + 1;
            nYGoBack = nY;
            nGoBackNewDirection = EAST;

            break;

         case NORTH_WEST:
            // The sea is towards the LHS (SW) of the coast, so first try to go left (to the SW)
            nXSeaward = nX - 1;
            nYSeaward = nY + 1;
            nSeawardNewDirection = SOUTH_WEST;

            // If can't do this, try to go straight on (to the NW)
            nXStraightOn = nX - 1;
            nYStraightOn = nY - 1;

            // If can't do either of these, try to go anti-seaward i.e. towards the RHS (NE)
            nXAntiSeaward = nX + 1;
            nYAntiSeaward = nY - 1;
            nAntiSeawardNewDirection = NORTH_EAST;

            // As a last resort, go back (to the SE)
            nXGoBack = nX + 1;
            nYGoBack = nY + 1;
            nGoBackNewDirection = SOUTH_EAST;

            break;
         }

         break;
      }

      // First, find out if any adjacent cell is flagged as sea
      bool const bOnCoast = nAdjacentCellIsSea(nX, nY);
      if (bOnCoast)
      {
         // Has the current cell already marked been marked as a coast cell (belonging to another coast)?
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
         {
            // Yes it is already marked as a coast cell
            // LogStream << m_ulIter << ":\t cell [" << nX << "][" << nY << "] already marked as coastline " << m_pRasterGrid->m_Cell[nX][nY].nGetCoastline() << endl;
         }
         else
         {
            // Not already marked, is this an intervention cell with the top above SWL?
            if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() >= m_dThisIterSWL))
            {
               // It is, so add it to the vector
               ILTempGridCRS.AppendIfNotPrevious(&Pti);
            }
            else if (m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() >= m_dThisIterSWL)
            {
               // The sediment top (inc any talus) is above SWL so add it to the vector object
               ILTempGridCRS.AppendIfNotPrevious(&Pti);
            }
         }
      }

      // Now find the next cell that we will move to: first try going in the direction of the sea. Is this seaward cell still within the grid?
      if (bIsWithinValidGrid(nXSeaward, nYSeaward))
      {
         // It is, so check if the cell in the seaward direction is a sea cell
         if (! m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsInContiguousSea())
         {
            // The seaward cell is not a sea cell, so we will move to it next time
            nX = nXSeaward;
            nY = nYSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nSeawardNewDirection;
            continue;
         }
      }

      // OK, we couldn't move seaward (but we may have marked the current cell as coast) so next try to move straight on. Is this straight-ahead cell still within the grid?
      if (bIsWithinValidGrid(nXStraightOn, nYStraightOn))
      {
         // It is, so check if there is sea immediately in front
         if (! m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsInContiguousSea())
         {
            // The straight-ahead cell is not a sea cell, so we will move to it next time
            nX = nXStraightOn;
            nY = nYStraightOn;

            // The search direction remains unchanged
            continue;
         }
      }

      // Couldn't move either seaward or straight on so next try to move in the anti-seaward direction. Is this anti-seaward cell still within the grid?
      if (bIsWithinValidGrid(nXAntiSeaward, nYAntiSeaward))
      {
         // It is, so check if there is sea at this anti-seaward cell
         if (! m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsInContiguousSea())
         {
            // The anti-seaward cell is not a sea cell, so we will move to it next time
            nX = nXAntiSeaward;
            nY = nYAntiSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nAntiSeawardNewDirection;
            continue;
         }
      }

      // Could not move to the seaward side, move straight ahead, or move to the anti-seaward side, so we must be in a single-cell dead end! As a last resort, turn round and move back to where we just came from, but first check that this is a valid cell
      if (bIsWithinValidGrid(nXGoBack, nYGoBack))
      {
         nX = nXGoBack;
         nY = nYGoBack;

         // And change the search direction
         nSearchDirection = nGoBackNewDirection;
      }
      else
      {
         // Our final choice is not a valid cell, so give up
         bDeadEnd = true;
         break;
      }

      // IF bOnCoast then MARK AS COAST

   } while (true);

   // OK, we have a possible coastline
   int nCoastSize = ILTempGridCRS.nGetSize();

   // Safety check
   if (nCoastSize == 0)
   {
      // Zero length vector coastline, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ":\t abandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since is zero length" << endl;

      return RTN_ERR_COAST_TRACING_ZERO_LENGTH;
   }

   // OK so far
   int nEndX = ILTempGridCRS[nCoastSize - 1].nGetX();
   int nEndY = ILTempGridCRS[nCoastSize - 1].nGetY();

   // Do the first check
   if ((nStartX == nEndX) && (nStartY == nEndY))
   {
      // Coastline starts and ends at same cell, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": abandoning possible coastline since it both starts and finishes at [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;

      return RTN_ERR_COAST_TRACING_SAME_START_FINISH;
   }

   // Do other checks
   if (bDeadEnd)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ":\t abandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since hit dead end at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, coastline size is " << nCoastSize << endl;

      return RTN_ERR_COAST_TRACING_OFF_EDGE;
   }

   // if (bTooLong)
   // {
   //    // Around loop too many times, so abandon this coastline
   //    if (m_nLogFileDetail >= LOG_FILE_ALL)
   //    {
   //       LogStream << m_ulIter << ":\t abandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since round loop " << nRoundLoop << " times, coastline size is " << nCoastSize;
   //
   //       if (nCoastSize > 0)
   //          LogStream << ", ended at [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "}";
   //
   //       LogStream << endl;
   //    }
   //
   //    return RTN_ERR_COAST_TRACING_TOO_LONG;
   // }

   if (bHereBefore)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         LogStream << m_ulIter << ":\t abandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since have visited cel before, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "}";

         LogStream << endl;
      }

      return RTN_ERR_COAST_TRACING_REPEATING;
   }

   if (nCoastSize < m_nCoastMin)
   {
      // The vector coastline is too small, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ":\t abandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;

      return RTN_ERR_COAST_TRACING_TOO_SHORT;
   }

   // Check the last point
   if ((nEndX != nX) || (nEndY != nY))
   {
      // Is the final cell in ILTempGridCRS already at the edge of the grid?
      if (! m_pRasterGrid->m_Cell[nEndX][nEndY].bIsBoundingBoxEdge())
      {
         // The final cell in ILTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
         ILTempGridCRS.AppendIfNotPrevious(nX, nY);
         nCoastSize++;

         nEndX = nX;
         nEndY = nY;
      }
   }

   // Get the grid edge at which the coastline finishes
   int const nEndEdge = m_pRasterGrid->m_Cell[nEndX][nEndY].nGetBoundingBoxEdge();

   if (nStartEdge == nEndEdge)
   {
      // Coastline starts and ends at same edge of grid, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": abandoning possible coastline since it both starts and finishes at [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;

      return RTN_ERR_COAST_TRACING_SAME_EDGE_START_FINISH;
   }

   // We have a valid coastline. Flag the coast start and end cells so that they are no longer possible coast start cells
   for (int nn = 0; nn < static_cast<int>(pV2DIPossibleStartCell->size()); nn++)
   {
      if (nn == nTraceFromStartCellIndex)
      {
         // Current start cell
         pVbTraced->at(nn) = true;
         continue;
      }

      int const nPossStartCellX = pV2DIPossibleStartCell->at(nn).nGetX();
      int const nPossStartCellY = pV2DIPossibleStartCell->at(nn).nGetY();

      // Need to check cells on either side, since coast finish cell can be offset by one cell
      int nXTmp;
      int nYTmp;
      for (int mmX = -1; mmX <= 1; mmX++)
      {
         for (int mmY = -1; mmY <= 1; mmY++)
         {
            nXTmp = nEndX + mmX;
            nYTmp = nEndY + mmY;

            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if ((nXTmp == nPossStartCellX) && (nYTmp == nPossStartCellY))
                  // Current end cell
                  pVbTraced->at(nn) = true;
            }
         }
      }
   }

   // Next, convert the grid coordinates in ILTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
   CGeomLine LTempExtCRS;

   for (int j = 0; j < nCoastSize; j++)
      LTempExtCRS.Append(dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()));

   // Now do some smoothing of the vector output, if desired
   if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
      LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS);
   else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);
   else if (m_nCoastSmooth == SMOOTH_RUNNING_MEDIAN)
      LTempExtCRS = LSmoothCoastRunningMedian(&LTempExtCRS);

   //    // DEBUG CODE ==================================================================================================
   // LogStream << "==================================" << endl;
   // for (int j = 0; j < nCoastSize; j++)
   // {
   // LogStream << "{" << dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()) << "}" << "\t{" << LTempExtCRS.dGetXAt(j) << ", " << LTempExtCRS.dGetYAt(j) << "}" << endl;
   // }
   // LogStream << "==================================" << endl;
   //    // DEBUG CODE ==================================================================================================

   // Create a new coastline object and append to it the vector of coastline objects
   CRWCoast const CoastTmp(this);
   m_VCoast.push_back(CoastTmp);
   int const nCoast = static_cast<int>(m_VCoast.size()) - 1;

   // Now mark the coastline on the grid
   for (int n = 0; n < nCoastSize; n++)
      m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsCoastline(nCoast);

   // Set the coastline (Ext CRS)
   m_VCoast[nCoast].SetCoastlineExtCRS(&LTempExtCRS);

   // Set the coastline (Grid CRS)
   m_VCoast[nCoast].SetCoastlineGridCRS(&ILTempGridCRS);

   // CGeom2DPoint PtLast(DBL_MIN, DBL_MIN);
   // for (int j = 0; j < nCoastSize; j++)
   // {
   //       // Store the smoothed points (in external CRS) in the coast's m_LCoastlineExtCRS object, also append dummy values to the other attribute vectors
   // if (PtLast != &LTempExtCRS[j])        // Avoid duplicate points
   // {
   // m_VCoast[nCoast].AppendPointToCoastlineExtCRS(LTempExtCRS[j].dGetX(), LTempExtCRS[j].dGetY());
   //
   //          // Also store the locations of the corresponding unsmoothed points (in raster grid CRS) in the coast's m_ILCellsMarkedAsCoastline vector
   // m_VCoast[nCoast].AppendCellMarkedAsCoastline(&ILTempGridCRS[j]);
   // }
   //
   // PtLast = LTempExtCRS[j];
   // }

   // Set values for the coast's other attributes: set the coast's handedness, and start and end edges
   m_VCoast[nCoast].SetSeaHandedness(nHandedness);
   m_VCoast[nCoast].SetStartEdge(nStartEdge);
   m_VCoast[nCoast].SetEndEdge(nEndEdge);

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
   {
      LogStream << m_ulIter << ":\t valid coast " << nCoast << " created, from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "} with " << nCoastSize << " points, handedness = " << (nHandedness == LEFT_HANDED ? "left" : "right") << endl;

      LogStream << m_ulIter << ":\t smoothed coastline " << nCoast << " runs from {" << LTempExtCRS[0].dGetX() << ", " << LTempExtCRS[0].dGetY() << "} to {" << LTempExtCRS[nCoastSize - 1].dGetX() << ", " << LTempExtCRS[nCoastSize - 1].dGetY() << "} i.e. from the ";

      if (nStartEdge == NORTH)
         LogStream << "north";
      else if (nStartEdge == SOUTH)
         LogStream << "south";
      else if (nStartEdge == WEST)
         LogStream << "west";
      else if (nStartEdge == EAST)
         LogStream << "east";

      LogStream << " edge to the ";
      if (nEndEdge == NORTH)
         LogStream << "north";
      else if (nEndEdge == SOUTH)
         LogStream << "south";
      else if (nEndEdge == WEST)
         LogStream << "west";
      else if (nEndEdge == EAST)
         LogStream << "east";
      LogStream << " edge" << endl;
   }

   // LogStream << "-----------------" << endl;
   // for (int kk = 0; kk < m_VCoast.back().nGetCoastlineSize(); kk++)
   // LogStream << kk << " [" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX() << "][" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY()) << "}" << endl;
   // LogStream << "-----------------" << endl;

   // Next calculate the curvature of the vector coastline
   DoCoastCurvature(nCoast, nHandedness);

   // Calculate values for the coast's flux orientation vector
   CalcCoastTangents(nCoast);

   // And create the vector of pointers to coastline-normal objects
   m_VCoast[nCoast].CreateProfilesAtCoastPoints();

   return RTN_OK;}

//===============================================================================================================================
//! Returns true if any of the eight surrounding cells is flagged as sea
//===============================================================================================================================
bool CSimulation::nAdjacentCellIsSea(int const nX, int const nY)
{
   int nXAdj;
   int nYAdj;

   for (int nSearchDirection = NORTH; nSearchDirection <= NORTH_WEST; nSearchDirection++)
   {
      switch (nSearchDirection)
      {
      case NORTH:
         nXAdj = nX;
         nYAdj = nY - 1;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case NORTH_EAST:
         nXAdj = nX + 1;
         nYAdj = nY - 1;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case EAST:
         // The sea is towards the RHS (S) of the coast, so first try to go right (to the S)
         nXAdj = nX + 1;
         nYAdj = nY;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case SOUTH_EAST:
         nXAdj = nX + 1;
         nYAdj = nY + 1;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case SOUTH:
         // The sea is towards the RHS (W) of the coast, so first try to go right (to the W)
         nXAdj = nX;
         nYAdj = nY + 1;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case SOUTH_WEST:
         nXAdj = nX - 1;
         nYAdj = nY + 1;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case WEST:
         // The sea is towards the RHS (N) of the coast, so first try to go right (to the N)
         nXAdj = nX - 1;
         nYAdj = nY;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;

      case NORTH_WEST:
         nXAdj = nX - 1;
         nYAdj = nY - 1;
         if (bIsWithinValidGrid(nXAdj, nYAdj))
         {
            if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
               return true;
         }
         break;
      }
   }

   return false;
}
