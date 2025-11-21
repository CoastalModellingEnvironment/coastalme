/*!
   \file locate_coast.cpp
   \brief Finds the coastline on the raster grid
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
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

#include <cstdlib>
using std::exit;

#include <string>
using std::to_string;

#include <gdal.h>
#include <gdal_priv.h>
#include <cpl_string.h>

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include <ios>
using std::fixed;

#include <iomanip>
using std::setprecision;

#include <stack>
using std::stack;

#include "cme.h"
#include "2di_point.h"
#include "i_line.h"
#include "line.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! First find all inundated cells and coast cells, then locate the vector coastline(s)
//===============================================================================================================================
int CSimulation::nLocateSeaAndCoasts(void)
{
   // Find all connected sea cells
   FindAllSeaCellsAndMarkCoastCells();

   // Create the vector coastline(s)
   int const nRet = nTraceAllVectorCoasts();
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
//! Finds and flags all inundated cells, also flags all coast with the DUMMY_COAST value
//===============================================================================================================================
void CSimulation::FindAllSeaCellsAndMarkCoastCells(void)
{
   bool bIgnoreFirst = true;        // Flag used to ignore the first X value in each X row
   bool bLastCellIsSea = false;
   bool bThisCellIsSea = false;

   // Go through all cells: X rows first, then Y columns
   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
         {
            m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();

            // Set this sea cell to have deep water (off-shore) wave orientation and height, will change this later for cells closer to the shoreline if we have on-shore waves
            m_pRasterGrid->m_Cell[nX][nY].SetWaveValuesToDeepWaterWaveValues();

            bThisCellIsSea = true;
         }
         else
            bThisCellIsSea = false;

         if (bIgnoreFirst)
         {
            // We are at the start of the X row, so ignore the last "Is it sea?" value
            bIgnoreFirst = false;
            bLastCellIsSea = bThisCellIsSea;
            continue;
         }

         // We are at the end of the X row, so set the switch ready for the next X value read
         if (nX == m_nXGridSize-1)
            bIgnoreFirst = true;

         if (bThisCellIsSea && (! bLastCellIsSea))
         {
            // This cell is sea and the previous cell was not sea, so flag the previous cell to be coast
            m_pRasterGrid->m_Cell[nX-1][nY].SetAsCoastline(DUMMY_COAST);
         }

         if ((! bThisCellIsSea) && bLastCellIsSea)
         {
            // This cell is not sea, the previous cell was sea: so flag this cell as coast
            m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(DUMMY_COAST);
         }

         bLastCellIsSea = bThisCellIsSea;
      }
   }

   bIgnoreFirst = true;        // Now used to ignore the first Y value in each Y column
   bLastCellIsSea = false;
   bThisCellIsSea = false;

   // Again: go through all cells: this time, Y columns first then X rows
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
            bThisCellIsSea = true;
         else
            bThisCellIsSea = false;

         if (bIgnoreFirst)
         {
            // We are at the start of the Y column, so ignore the last "Is it sea?" value
            bIgnoreFirst = false;
            bLastCellIsSea = bThisCellIsSea;
            continue;
         }

         // We are at the end of the Y column, so set the switch ready for the next Y value read
         if (nY == m_nYGridSize-1)
            bIgnoreFirst = true;

         if (bThisCellIsSea && (! bLastCellIsSea))
         {
            // This cell is sea and the previous cell was not sea, so flag the previous cell to be coast
            m_pRasterGrid->m_Cell[nX][nY-1].SetAsCoastline(DUMMY_COAST);
         }

         if ((! bThisCellIsSea) && bLastCellIsSea)
         {
            // This cell is not sea, the previous cell was sea: so flag this cell as coast
            m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(DUMMY_COAST);
         }

         bLastCellIsSea = bThisCellIsSea;
      }
   }

   // // DEBUG CODE ===========================================================================================================
   // // string strOutFile1 = m_strOutPath + "is_contiguous_sea_";
   // // string strOutFile1 = m_strOutPath + "is_inundated_";
   // string strOutFile1 = m_strOutPath + "is_coastline_";
   // strOutFile1 += to_string(m_ulIter);
   // strOutFile1 += ".tif";
   //
   // GDALDriver* pDriver1 = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet1 = pDriver1->Create(strOutFile1.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet1->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet1->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster1 = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   //    for (int nX = 0; nX < m_nXGridSize; nX++)
   //    {
   //       // pdRaster1[n++] = m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea();
   //       // pdRaster1[n++] = m_pRasterGrid->m_Cell[nX][nY].bIsInundated();
   //       pdRaster1[n++] = m_pRasterGrid->m_Cell[nX][nY].nGetCoastline();
   //    }
   // }
   //
   // GDALRasterBand* pBand1 = pDataSet1->GetRasterBand(1);
   // pBand1->SetNoDataValue(m_dMissingValue);
   // int nRet1 = pBand1->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster1, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet1 == CE_Failure)
   //    exit(1);
   //
   // GDALClose(pDataSet1);
   // delete[] pdRaster1;
   // // DEBUG CODE ===========================================================================================================
}

//===============================================================================================================================
//! Locates all coastline start/finish points on the edges of the raster grid, then traces vector coastline(s) from these start points
//===============================================================================================================================
int CSimulation::nTraceAllVectorCoasts(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Tracing vector coasts" << endl;

   int nValidCoast = 0;
   int nNumEdgeCell = static_cast<int>(m_VEdgeCell.size());
   vector<bool> VbTraced(nNumEdgeCell, false);

   // Go along the list of edge cells and look for coastline start/finish cells
   for (int nEdgeCell = 0; nEdgeCell < static_cast<int>(m_VEdgeCell.size()); nEdgeCell++)
   {
      if (m_bOmitSearchNorthEdge && (m_VEdgeCellEdge[nEdgeCell] == NORTH || m_VEdgeCellEdge[nEdgeCell + 1] == NORTH))
         continue;

      if (m_bOmitSearchSouthEdge && (m_VEdgeCellEdge[nEdgeCell] == SOUTH || m_VEdgeCellEdge[nEdgeCell + 1] == SOUTH))
         continue;

      if (m_bOmitSearchWestEdge && (m_VEdgeCellEdge[nEdgeCell] == WEST || m_VEdgeCellEdge[nEdgeCell + 1] == WEST))
         continue;

      if (m_bOmitSearchEastEdge && (m_VEdgeCellEdge[nEdgeCell] == EAST || m_VEdgeCellEdge[nEdgeCell + 1] == EAST))
         continue;

      int const nX = m_VEdgeCell[nEdgeCell].nGetX();
      int const nY = m_VEdgeCell[nEdgeCell].nGetY();

      if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline() && (! VbTraced[nEdgeCell]))
      {
         // This edge cell could be the start of a new coastline
         int nHandedness;

         // Calculate handedness of this potential coastline (shows which side the sea is on when travelling down-coast i.e. in the direction in which coastline point numbers INCREASE)
         int nStartEdge = m_VEdgeCellEdge[nEdgeCell];
         if (nStartEdge == NORTH)
         {
            if (nX > 0)
            {
               if (m_pRasterGrid->m_Cell[nX-1][nY].bIsInContiguousSea())
                  nHandedness = RIGHT_HANDED;
               else
                  nHandedness = LEFT_HANDED;
            }
            else if (nX < m_nXGridSize-1)
            {
               if (m_pRasterGrid->m_Cell[nX+1][nY].bIsInContiguousSea())
                  nHandedness = LEFT_HANDED;
               else
                  nHandedness = RIGHT_HANDED;
            }
         }
         else if (nStartEdge == SOUTH)
         {
            if (nX > 0)
            {
               if (m_pRasterGrid->m_Cell[nX-1][nY].bIsInContiguousSea())
                  nHandedness = LEFT_HANDED;
               else
                  nHandedness = RIGHT_HANDED;
            }
            else if (nX < m_nXGridSize-1)
            {
               if (m_pRasterGrid->m_Cell[nX+1][nY].bIsInContiguousSea())
                  nHandedness = RIGHT_HANDED;
               else
                  nHandedness = LEFT_HANDED;
            }
         }
         else if (nStartEdge == WEST)
         {
            if (nY > 0)
            {
               if (m_pRasterGrid->m_Cell[nX][nY-1].bIsInContiguousSea())
                  nHandedness = LEFT_HANDED;
               else
                  nHandedness = RIGHT_HANDED;
            }
            else if (nY < m_nYGridSize-1)
            {
               if (m_pRasterGrid->m_Cell[nX][nY+1].bIsInContiguousSea())
                  nHandedness = RIGHT_HANDED;
               else
                  nHandedness = LEFT_HANDED;
            }
         }
         else if (nStartEdge == EAST)
         {
            if (nY > 0)
            {
               if (m_pRasterGrid->m_Cell[nX][nY-1].bIsInContiguousSea())
                  nHandedness = RIGHT_HANDED;
               else
                  nHandedness = LEFT_HANDED;
            }
            else if (nY < m_nYGridSize-1)
            {
               if (m_pRasterGrid->m_Cell[nX][nY+1].bIsInContiguousSea())
                  nHandedness = LEFT_HANDED;
               else
                  nHandedness = RIGHT_HANDED;
            }
         }

         CGeom2DIPoint PtiStart(nX, nY);
         CGeom2DIPoint PtiEnd;

         int nRet = nTraceVectorCoastLine(nStartEdge, nHandedness, &PtiStart, &PtiEnd);
         if (nRet == RTN_OK)
         {
            // We have a valid coastline starting from this possible start cell
            VbTraced[nEdgeCell] = true;

            // Find the end cell in the list of edge cells
            auto it = find(m_VEdgeCell.begin(), m_VEdgeCell.end(), &PtiEnd);

            if (it == m_VEdgeCell.end())
            {
               // Not found. This can happen because of rounding problems, i.e. the cell which was stored as the first cell of the raster coastline
               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  LogStream << m_ulIter << ": " << ERR << " when tracing coast, [" << PtiEnd.nGetX() << "][" << PtiEnd.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiEnd.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiEnd.nGetY()) << "} not found in list of edge cells" << endl;

               return RTN_ERR_COAST_CANT_FIND_EDGE_CELL;
            }

            // Found
            int nPos = static_cast<int>(it - m_VEdgeCell.begin());
            VbTraced[nPos] = true;

            nValidCoast++;
         }
      }
   }

   // Any possible coastline start/finish cells found?
   if (nValidCoast == 0)
   {
      LogStream << m_ulIter << ": no coastline start/finish points found after grid edges searched.";

      if (m_bOmitSearchNorthEdge || m_bOmitSearchSouthEdge || m_bOmitSearchWestEdge || m_bOmitSearchEastEdge)
      {
         LogStream << " Note that the following grid edges were not searched: " << (m_bOmitSearchNorthEdge ? "N " : "") << (m_bOmitSearchSouthEdge ? "S " : "") << (m_bOmitSearchWestEdge ? "W " : "") << (m_bOmitSearchEastEdge ? "E " : "");
      }

      LogStream << endl;

      return RTN_ERR_NO_START_FINISH_POINTS_TRACING_COAST;
   }

   // if (nValidCoast == 0)
   // {
   //    // No valid coasts found so try again, this time working through the possible start/finish points in reverse order
   //    for (int n = 0; n < static_cast<int>(VbTraced.size()); n++)
   //       VbTraced[n] = false;
   //
   //    for (int n = 0; n < static_cast<int>(V2DIPossibleStartCell.size()); n++)
   //    {
   //       if (! VbTraced[n])
   //       {
   //          int nRet = 0;
   //
   //          if (VbPossibleStartCellLHEdge[n])
   //             nRet = nTraceVectorCoastLine(LEFT_HANDED, &V2DIPossibleStartCell);
   //          else
   //             nRet = nTraceVectorCoastLine(RIGHT_HANDED, &V2DIPossibleStartCell);
   //
   //          if (nRet == RTN_OK)
   //          {
   //             // We have a valid coastline starting from this possible start cell
   //             VbTraced[n] = true;
   //             nValidCoast++;
   //          }
   //       }
   //    }
   // }
   //
   if (nValidCoast == 0)
   {
      // Still no valid coasts found, so we have to give up
      cerr << m_ulIter << ": no valid coasts found, see " << m_strLogFile << " for more information" << endl;
      return RTN_ERR_NO_VALID_COAST;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Traces a coastline which is marked on the raster grid from the given grid edge to another grid edge. The resulting vector coastline is then smoothed TODO Only from edges, need to deal with inland water
//===============================================================================================================================
int CSimulation::nTraceVectorCoastLine(int const nStartEdge, int const nHandedness, CGeom2DIPoint const* pPtiStartCell, CGeom2DIPoint* pPtiEndCell)
{
   LogStream << m_ulIter << ": \ttracing coastline from ";
   switch(nStartEdge)
   {
      case NORTH:
         LogStream << "NORTH";
         break;
      case NORTH_EAST:
         LogStream << "NORTH EAST";
         break;
      case EAST:
         LogStream << "EAST";
         break;
      case SOUTH_EAST:
         LogStream << "SOUTH EAST";
         break;
      case SOUTH:
         LogStream << "SOUTH";
         break;
      case SOUTH_WEST:
         LogStream << "SOUTH WEST";
         break;
      case WEST:
         LogStream << "WEST";
         break;
      case NORTH_WEST:
         LogStream << "NORTH WEST";
         break;
   }
   LogStream << " edge, handedness = ";
   if (nHandedness == LEFT_HANDED)
      LogStream << "LEFT-HANDED";
   else
      LogStream << "RIGHT HANDED";
   LogStream << endl;

   // Temporary coastline as integer points (grid CRS)
   CGeomILine ILTempGridCRS;

   // Add the start cell to the vector
   ILTempGridCRS.Append(pPtiStartCell);

   bool bHitEdge = false;
   bool bTooLong = false;
   bool bRepeating = false;

   int nRoundLoop = -1;
   int nX = pPtiStartCell->nGetX();
   int nY = pPtiStartCell->nGetY();
   int nStartX = nX;
   int nStartY = nY;

   do
   {
      bool bFound = false;

      for (int nSearchDirection = NORTH; nSearchDirection <= NORTH_WEST; nSearchDirection++)
      {
         if (bFound)
            break;

         int nXAdj;
         int nYAdj;

         switch (nSearchDirection)
         {
         case NORTH:
            nXAdj = nX;
            nYAdj = nY - 1;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to NORTH marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case NORTH_EAST:
            nXAdj = nX + 1;
            nYAdj = nY - 1;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to NORTH EAST marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case EAST:
            nXAdj = nX + 1;
            nYAdj = nY;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to EAST marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case SOUTH_EAST:
            nXAdj = nX + 1;
            nYAdj = nY + 1;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to SOUTH EAST marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case SOUTH:
            nXAdj = nX;
            nYAdj = nY + 1;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to SOUTH marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case SOUTH_WEST:
            nXAdj = nX - 1;
            nYAdj = nY + 1;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to SOUTH WEST marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case WEST:
            nXAdj = nX - 1;
            nYAdj = nY;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tcell to WEST marked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;

         case NORTH_WEST:
            nXAdj = nX - 1;
            nYAdj = nY - 1;

            if (bIsWithinValidGrid(nXAdj, nYAdj))
            {
               CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
               if (! ILTempGridCRS.bIsPresent(&PtiTmp))
               {
                  int nCoastline = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline();
                  if (nCoastline == DUMMY_COAST)
                  {
                     // This cell was marked as a coastline
                     ILTempGridCRS.Append(&PtiTmp);
                     nX = nXAdj;
                     nY = nYAdj;
                     bFound = true;
                     LogStream << m_ulIter << ": \tmarked as coastline [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                     break;
                  }
               }
            }
            break;
         }
      }

      // Did we hit an edge?
      if (m_pRasterGrid->m_Cell[nX][nY].bIsBoundingBoxEdge())
      {
          bHitEdge = true;
          LogStream << m_ulIter << ": \thit edge at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

          break;
      }

      if (! bFound)
      {
         LogStream << m_ulIter << ": \tnot found at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         break;
      }

      // Safety check
      if (++nRoundLoop > m_nCoastMax)
      {
         bTooLong = true;
         LogStream << m_ulIter << ": \ttoo long [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         break;
      }

      // Another safety check
      if ((nRoundLoop > 10) && (ILTempGridCRS.nGetSize() < 2))
      {
         // We've been 10 times round the loop but the coast is still less than 2 coastline points in length, so we must be repeating
         bRepeating = true;
         LogStream << m_ulIter << ": \trepeating at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

         break;
      }

   } while (true);

   // OK, we have a coastline. But is it any good?
   int nCoastSize = ILTempGridCRS.nGetSize();

   // Now check this possible coastline
   if (bHitEdge && (*pPtiStartCell == ILTempGridCRS.Back()))
   {
      // This coastline is a loop, returning to the start cell
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      {
         LogStream << m_ulIter << ": \tabandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since returned to start point after looping " << nRoundLoop << " times, coastline size is " << nCoastSize << endl;
      }

      // Unmark the cells
      for (int nCell = 0; nCell < nCoastSize; nCell++)
         m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(INT_NODATA);

      return RTN_ERR_IGNORING_COAST;
   }

   if (bTooLong)
   {
      // Around loop too many times, so abandon this coastline
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      {
         LogStream << m_ulIter << ": \tabandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since round loop " << nRoundLoop << " times, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      // Unmark the cells
      for (int nCell = 0; nCell < nCoastSize; nCell++)
         m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(INT_NODATA);

      return RTN_ERR_TOO_LONG_TRACING_COAST;
   }

   if (bRepeating)
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      {
         LogStream << m_ulIter << ": abandon possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since repeating, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      // Unmark the cells
      for (int nCell = 0; nCell < nCoastSize; nCell++)
         m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(INT_NODATA);

      return RTN_ERR_REPEATING_WHEN_TRACING_COAST;
   }

   if (nCoastSize < m_nCoastMin)
   {
      // The vector coastline is too small, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": \tabandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;

      // Unmark the cells
      for (int nCell = 0; nCell < nCoastSize; nCell++)
         m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(INT_NODATA);

      return RTN_ERR_COAST_TOO_SMALL;
   }

   // OK this new coastline is fine
   int const nEndX = nX;
   int const nEndY = nY;
   int const nCoastEndX = ILTempGridCRS[nCoastSize - 1].nGetX();
   int const nCoastEndY = ILTempGridCRS[nCoastSize - 1].nGetY();

   if ((nCoastEndX != nEndX) || (nCoastEndY != nEndY))
   {
      // The grid-edge cell at nEndX, nEndY is not already at end of ILTempGridCRS. But is the final cell in ILTempGridCRS already at the edge of the grid?
      if (! m_pRasterGrid->m_Cell[nCoastEndX][nCoastEndY].bIsBoundingBoxEdge())
      {
         // The final cell in ILTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
         ILTempGridCRS.AppendIfNotPrevious(nEndX, nEndY);
         nCoastSize++;
      }
   }

   pPtiEndCell->SetX(nEndX);
   pPtiEndCell->SetY(nEndY);

   // Need to also specify end edge for smoothing routines
   int const nEndEdge = m_pRasterGrid->m_Cell[nEndX][nEndY].nGetBoundingBoxEdge();

   // Next, convert the grid coordinates in ILTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
   CGeomLine LTempExtCRS;

   for (int j = 0; j < nCoastSize; j++)
      LTempExtCRS.Append(dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()));

   // Now do some smoothing of the vector output, if desired
   if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
      LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS);
   else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);

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
      LogStream << m_ulIter << ": \tvalid coast " << nCoast << " created, from [" << nStartX << "][" << nStartY << "] to [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "} with " << nCoastSize << " points, handedness = " << (nHandedness == LEFT_HANDED ? "left" : "right") << endl;

      LogStream << m_ulIter << ": \tsmoothed coastline " << nCoast << " runs from {" << LTempExtCRS[0].dGetX() << ", " << LTempExtCRS[0].dGetY() << "} to {" << LTempExtCRS[nCoastSize - 1].dGetX() << ", " << LTempExtCRS[nCoastSize - 1].dGetY() << "} i.e. from the ";

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

   return RTN_OK;
}

//===============================================================================================================================
//! First find all connected sea areas, then locate the vector coastline(s), then put these onto the raster grid
//===============================================================================================================================
int CSimulation::nLocateFloodAndCoasts(void)
{
   // Find all connected sea cells
   FindAllInundatedCells();

   // Find every coastline on the raster grid, mark raster cells, then create the vector coastline
   int const nRet = nTraceAllFloodCoasts();

   if (nRet != RTN_OK)
      return nRet;

   // Have we created any coasts?
   switch (m_nLevel)
   {
   case 0: // WAVESETUP + SURGE:
   {
      if (m_VFloodWaveSetupSurge.empty())
      {
         cerr << m_ulIter << ": " << ERR << "no flood coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
         return RTN_ERR_NO_COAST;
      }

      break;
   }

   case 1: // WAVESETUP + SURGE + RUNUP:
   {
      if (m_VFloodWaveSetupSurgeRunup.empty())
      {
         cerr << m_ulIter << ": " << ERR << "no flood coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
         return RTN_ERR_NO_COAST;
      }

      break;
   }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Finds and flags all sea areas which have at least one cell at a grid edge (i.e. does not flag 'inland' seas)
//===============================================================================================================================
int CSimulation::FindAllInundatedCells(void)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         m_pRasterGrid->m_Cell[nX][nY].UnSetCheckFloodCell();
         m_pRasterGrid->m_Cell[nX][nY].UnSetInContiguousFlood();
         m_pRasterGrid->m_Cell[nX][nY].SetAsFloodline(false);
      }
   }

   // Go along the list of edge cells
   for (unsigned int n = 0; n < m_VEdgeCell.size(); n++)
   {
      if (m_bOmitSearchNorthEdge && m_VEdgeCellEdge[n] == NORTH)
         continue;

      if (m_bOmitSearchSouthEdge && m_VEdgeCellEdge[n] == SOUTH)
         continue;

      if (m_bOmitSearchWestEdge && m_VEdgeCellEdge[n] == WEST)
         continue;

      if (m_bOmitSearchEastEdge && m_VEdgeCellEdge[n] == EAST)
         continue;

      int const nX = m_VEdgeCell[n].nGetX();
      int const nY = m_VEdgeCell[n].nGetY();

      if ((! m_pRasterGrid->m_Cell[nX][nY].bIsCellFloodCheck()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()))
      {
         // This edge cell is below SWL and sea depth remains set to zero
         FloodFillLand(nX, nY);
      }
   }

   return RTN_OK;
}
