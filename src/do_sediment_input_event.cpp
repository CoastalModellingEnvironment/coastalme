/*!
   \file do_sediment_input_event.cpp
   \brief Deposits sediment onto the grid
   \details TODO 001 A more detailed description of these routines.
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
#include <cstdio>

#include <iostream>
using std::endl;

#include <algorithm>
using std::begin;
using std::end;
using std::find;

#include "cme.h"
#include "simulation.h"
#include "cell.h"
#include "coast.h"
#include "sediment_input_event.h"
#include "2di_point.h"

//===============================================================================================================================
//! Check to see if we have any sediment input events this timestep, if so then do the event(s)
//===============================================================================================================================
int CSimulation::nCheckForSedimentInputEvent(void)
{
   m_bSedimentInputThisIter = false;

   // Go through all sediment input events, check for any this timestep
   int const nEvents = static_cast<int>(m_pVSedInputEvent.size());

   for (int n = 0; n < nEvents; n++)
   {
      if (m_pVSedInputEvent[n]->ulGetEventTimeStep() == m_ulIter)
      {
         m_bSedimentInputThisIter = true;

         int const nRet = nDoSedimentInputEvent(n);

         if (nRet != RTN_OK)
            return nRet;
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Do a sediment input event
//===============================================================================================================================
int CSimulation::nDoSedimentInputEvent(int const nEvent)
{
   // Get values for the sediment input event
   int const nLocID = m_pVSedInputEvent[nEvent]->nGetLocationID();
   double const dFineSedVol = m_pVSedInputEvent[nEvent]->dGetFineSedVol();
   double const dSandSedVol = m_pVSedInputEvent[nEvent]->dGetSandSedVol();
   double const dCoarseSedVol = m_pVSedInputEvent[nEvent]->dGetCoarseSedVol();
   double const dLen = m_pVSedInputEvent[nEvent]->dGetLen();
   double const dWidth = m_pVSedInputEvent[nEvent]->dGetWidth();
   // double dThick = m_pVSedInputEvent[nEvent]->dGetThick();

   // TODO 083 Get all three kinds of sediment input events working correctly

   if (m_bSedimentInputAtPoint || m_bSedimentInputAtCoast)
   {
      // The sediment input event is at a user-specified point, or in a block at the nearest point on a coast to a user-specified point. So get the point from values read from the shapefile
      int const nEvents = static_cast<int>(m_VnSedimentInputLocationID.size());
      int nPointGridX = -1;
      int nPointGridY = -1;

      for (int n = 0; n < nEvents; n++)
      {
         if (m_VnSedimentInputLocationID[n] == nLocID)
         {
            nPointGridX = nRound(m_VdSedimentInputLocationX[n]); // In grid coords
            nPointGridY = nRound(m_VdSedimentInputLocationY[n]); // In grid coords
         }
      }

      if (nPointGridX == -1)
         // Should never get here
         return RTN_ERR_SEDIMENT_INPUT_EVENT;

      // All OK, so get landform
      CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLandform();

      // Is this sediment input event at a pre-specified fixed point, or in a block on a coast, or where a line intersects with a coast?
      if (m_bSedimentInputAtPoint)
      {
         // Sediment input is at a user-specified point
         int const nTopLayer = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].nGetNumOfTopLayerAboveBasement();

         // Is this user-specified point in a polygon?
         int const nThisPoly = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].nGetPolygonID();
         int nThisPolyCoast = INT_NODATA;
         if (nThisPoly != INT_NODATA)
         {
            // Yes we are in a polygon, so get the coast ID of the polygon for this cell
            nThisPolyCoast = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].nGetPolygonCoastID();

            // Safety check
            if (nThisPolyCoast == INT_NODATA)
               return RTN_ERR_CELL_IN_POLY_BUT_NO_POLY_COAST;
         }

         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         {
            LogStream << m_ulIter << ": Sediment input event " << nEvent + 1 << " at point [" << nPointGridX << "][" << nPointGridY << "] = {" << dGridXToExtCRSX(nPointGridX) << ", " << dGridYToExtCRSY(nPointGridY) << "] with location ID " << nLocID;

            if (nThisPoly != INT_NODATA)
               LogStream << " which is within coast " << nThisPolyCoast << " polygon " << nThisPoly << endl;
            else
               LogStream << " this is not within a polygon" << endl;
         }

         // Is some fine unconsolidated sediment being input?
         double const dFineDepth = dFineSedVol / m_dCellArea;
         if (dFineDepth > 0)
         {
            // Yes, so add to this cell's fine unconsolidated sediment
            m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFineSedimentInputDepth(dFineDepth);

            // And update the sediment top elevation value
            m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].CalcAllLayerElevsAndD50();

            // If we are in a polygon, then add to this polygon's sand sediment input total
            if (nThisPoly != INT_NODATA)
               m_VCoast[nThisPolyCoast].pGetPolygon(nThisPoly)->SetSedimentInputUnconsFine(dFineDepth);

            // Add to the this-iteration total of fine sediment input
            m_dThisiterUnconsFineInput += dFineDepth;

            // And assign the cell's landform category
            pLandform->SetLFCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
         }

         // Is some sand-sized unconsolidated sediment being input?
         double const dSandDepth = dSandSedVol / m_dCellArea;

         if (dSandDepth > 0)
         {
            // Yes, so add to this cell's sand unconsolidated sediment
            m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandSedimentInputDepth(dSandDepth);

            // And update the sediment top elevation value
            m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].CalcAllLayerElevsAndD50();

            // If we are in a polygon, then add to this polygon's sand sediment input total
            if (nThisPoly != INT_NODATA)
               m_VCoast[nThisPolyCoast].pGetPolygon(nThisPoly)->SetSedimentInputUnconsSand(dSandDepth);

            // Add to the this-iteration total of sand sediment input
            m_dThisiterUnconsSandInput += dSandDepth;

            // And assign the cell's landform category
            pLandform->SetLFCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
         }

         // Is some coarse unconsolidated sediment being input?
         double const dCoarseDepth = dCoarseSedVol / m_dCellArea;

         if (dCoarseDepth > 0)
         {
            // Yes, so add to this cell's coarse unconsolidated sediment
            m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseSedimentInputDepth(dCoarseDepth);

            // And update the sediment top elevation value
            m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].CalcAllLayerElevsAndD50();

            // If we are in a polygon, then add to this polygon's coarse sediment input total
            if (nThisPoly != INT_NODATA)
               m_VCoast[nThisPolyCoast].pGetPolygon(nThisPoly)->SetSedimentInputUnconsCoarse(dCoarseDepth);

            // Add to the this-iteration total of coarse sediment input
            m_dThisiterUnconsCoarseInput += dCoarseDepth;

            // And assign the cell's landform category
            pLandform->SetLFCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
         }

         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << ", depth of fine sediment added = " << dFineDepth << " m, depth of sand sediment added = " << dSandDepth << " m, depth of coarse sediment added = " << dCoarseDepth << " m" << endl;
      }
      else if (m_bSedimentInputAtCoast)
      {
         // Is in a sediment block, seaward from a coast
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": Sediment input event " << nEvent + 1 << " with location ID " << nLocID << " at closest point on coast to [" << nPointGridX << "][" << nPointGridY << "] = {" << dGridXToExtCRSX(nPointGridX) << ", " << dGridYToExtCRSY(nPointGridY) << "]" << endl;

         // Find the closest point on any coastline
         int nCoastClosest;
         CGeom2DIPoint const PtiCoastPoint = PtiFindClosestCoastPoint(nPointGridX, nPointGridY, nCoastClosest);

         int const nCoastX = PtiCoastPoint.nGetX();
         int const nCoastY = PtiCoastPoint.nGetY();

         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": Closest coast point is on coast " << nCoastClosest << " at [" << nCoastX << "][" << nCoastY << "] = {" << dGridXToExtCRSX(nCoastX) << ", " << dGridYToExtCRSY(nCoastY) << "}, along-coast width of sediment block = " << dWidth << " m, coast-normal length of sediment block = " << dLen << " m" << endl;

         int const nCoast = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetCoast();
         int const nCoastPoint = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetPointOnCoast();
         int const nHalfWidth = nRound(dWidth / m_dCellSide);
         int const nLength = nRound(dLen / m_dCellSide);
         int const nCoastLen = m_VCoast[nCoast].nGetCoastlineSize();
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

         int const nCoastHand = m_VCoast[nCoast].nGetSeaHandedness();
         int nPerpHand = LEFT_HANDED;

         if (nCoastHand == LEFT_HANDED)
            nPerpHand = RIGHT_HANDED;

         vector<CGeom2DIPoint> VPoints;

         // If size of plume is smaller than a cell, choose the current coast point
         if (nHalfWidth == 0)
         {
            VPoints.push_back(PtiCoastPoint);
         }

         else
         {
            // It is larger than a cell
            vector<int> VnCentrePointsXOffset, VnCentrePointsYOffset;

            for (int m = 0; m < nHalfWidth; m++)
            {
               if (m == 0)
               {
                  VPoints.push_back(CGeom2DIPoint(nCoastX, nCoastY));

                  for (int n = 1; n < nLength; n++)
                  {
                     CGeom2DIPoint const PtiTmp = PtiGetPerpendicular(nCoastXBefore, nCoastYBefore, nCoastXAfter, nCoastYAfter, n, nPerpHand);

                     if (bIsWithinValidGrid(&PtiTmp))
                     {
                        VPoints.push_back(PtiTmp);

                        VnCentrePointsXOffset.push_back(PtiTmp.nGetX() - nCoastX);
                        VnCentrePointsYOffset.push_back(PtiTmp.nGetY() - nCoastY);
                     }
                  }
               }
               else
               {
                  int const nCoastPointInBlockBefore = nCoastPoint - m;
                  int const nCoastPointInBlockAfter = nCoastPoint + m;

                  if (nCoastPointInBlockBefore >= 0)
                  {
                     int const nCoastXInBlockBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockBefore)->nGetX();
                     int const nCoastYInBlockBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockBefore)->nGetY();

                     if (bIsWithinValidGrid(nCoastXInBlockBefore, nCoastYInBlockBefore))
                     {
                        VPoints.push_back(CGeom2DIPoint(nCoastXInBlockBefore, nCoastYInBlockBefore));

                        for (unsigned int n = 0; n < VnCentrePointsXOffset.size(); n++)
                        {
                           int const nXTmp = nCoastXInBlockBefore + VnCentrePointsXOffset[n];
                           int const nYTmp = nCoastYInBlockBefore + VnCentrePointsYOffset[n];

                           if (bIsWithinValidGrid(nXTmp, nYTmp))
                              VPoints.push_back(CGeom2DIPoint(nXTmp, nYTmp));
                        }
                     }
                  }

                  if (nCoastPointInBlockAfter < nCoastLen)
                  {
                     int const nCoastXInBlockAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockAfter)->nGetX();
                     int const nCoastYInBlockAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockAfter)->nGetY();

                     if (bIsWithinValidGrid(nCoastXInBlockAfter, nCoastYInBlockAfter))
                     {
                        VPoints.push_back(CGeom2DIPoint(nCoastXInBlockAfter, nCoastYInBlockAfter));

                        for (unsigned int n = 0; n < VnCentrePointsXOffset.size(); n++)
                        {
                           int const nXTmp = nCoastXInBlockAfter + VnCentrePointsXOffset[n];
                           int const nYTmp = nCoastYInBlockAfter + VnCentrePointsYOffset[n];

                           if (bIsWithinValidGrid(nXTmp, nYTmp))
                           {
                              CGeom2DIPoint const PtiTmp(nXTmp, nYTmp);

                              // Is this cell already in the array?
                              if (find(begin(VPoints), end(VPoints), PtiTmp) == end(VPoints))
                                 // No it isn't
                                 VPoints.push_back(PtiTmp);
                           }
                        }
                     }
                  }
               }
            }

            // // DEBUG CODE ===============================================================================================================
            // LogStream << endl;
            // unsigned int m = 0;
            // for (unsigned int n = 0; n < VPoints.size(); n++)
            // {
            // LogStream << "[" << VPoints[n].nGetX() << ", " << VPoints[n].nGetY() << "] ";
            // m++;
            // if (m == VnCentrePointsXOffset.size() + 1)
            // {
            // LogStream << endl;
            // m = 0;
            // }
            // }
            // LogStream << endl;
            // // DEBUG CODE ===============================================================================================================
         }

         // OK we now know which cells are part of the sediment block, and so will receive sediment input. Next calculate the volume per cell
         double const dFineDepth = dFineSedVol / m_dCellArea;
         double const dSandDepth = dSandSedVol / m_dCellArea;
         double const dCoarseDepth = dCoarseSedVol / m_dCellArea;

         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": Total depth of fine sediment added = " << dFineDepth << " m, total depth of sand sediment added = " << dSandDepth << " m, total depth of coarse sediment added = " << dCoarseDepth << " m" << endl;

         size_t const nArea = VPoints.size();
         double const dArea = static_cast<double>(nArea);
         double const dFineDepthPerCell = dFineDepth / dArea;
         double const dSandDepthPerCell = dSandDepth / dArea;
         double const dCoarseDepthPerCell = dCoarseDepth / dArea;

         // OK, so finally: put some sediment onto each cell in the sediment block
         int const nTopLayer = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].nGetNumOfTopLayerAboveBasement();

         for (unsigned int n = 0; n < nArea; n++)
         {
            int const nX = VPoints[n].nGetX();
            int const nY = VPoints[n].nGetY();

            // Add to this cell's fine unconsolidated sediment
            m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFineSedimentInputDepth(dFineDepthPerCell);

            // Add to this cell's sand unconsolidated sediment
            m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandSedimentInputDepth(dSandDepthPerCell);

            // Add to this cell's coarse unconsolidated sediment
            m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseSedimentInputDepth(dCoarseDepthPerCell);
         }

         // Add to the this-iteration totals
         m_dThisiterUnconsFineInput += dFineDepth;
         m_dThisiterUnconsSandInput += dSandDepth;
         m_dThisiterUnconsCoarseInput += dCoarseDepth;
      }
   }
   else if (m_bSedimentInputAlongLine)
   {
      // The location of the sediment input event is where a line intersects a coast. First get the line, using values read from the shapefile
      int const nPoints = static_cast<int>(m_VnSedimentInputLocationID.size());
      vector<int> VnLineGridX;
      vector<int> VnLineGridY;

      for (int n = 0; n < nPoints; n++)
      {
         if (m_VnSedimentInputLocationID[n] == nLocID)
         {
            VnLineGridX.push_back(nRound(m_VdSedimentInputLocationX[n]));
            VnLineGridY.push_back(nRound(m_VdSedimentInputLocationY[n]));
         }
      }

      // // DEBUG CODE ===========================================================================================================
      // string strOutFile = m_strOutPath;
      // strOutFile += "00_sediment_input_line_CHECK_";
      // strOutFile += to_string(m_ulIter);
      // strOutFile += ".tif";
      //
      // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
      // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      // pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      // int nn = 0;
      // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
      // for (int nY = 0; nY < m_nYGridSize; nY++)
      // {
      // for (int nX = 0; nX < m_nXGridSize; nX++)
      // {
      // pdRaster[nn++] = 0;
      // }
      // }
      //
      // for (unsigned int n = 0; n < VnLineGridX.size() - 1; n++)
      // {
      // int nX = VnLineGridX[n];
      // int nY = VnLineGridY[n];
      // int m = (nY * m_nXGridSize) + nX;
      //
      // pdRaster[m] = 1;
      // }
      //
      // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
      // pBand->SetNoDataValue(m_dMissingValue);
      // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
      //
      // if (nRet == CE_Failure)
      // return RTN_ERR_GRIDCREATE;
      //
      // GDALClose(pDataSet);
      // delete[] pdRaster;
      // // DEBUG CODE ===========================================================================================================

      // Should never get here
      if (VnLineGridX.size() == 0)
         return RTN_ERR_SEDIMENT_INPUT_EVENT;

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ": Sediment input event " << nEvent << " at line/coast intersection for line with ID " << nLocID << endl;

      // Now go along the line, joining each pair of points by a straight line
      int nCoastX = -1;
      int nCoastY = -1;

      for (unsigned int n = 0; n < VnLineGridX.size() - 1; n++)
      {
         int const nXStart = VnLineGridX[n];
         int const nXEnd = VnLineGridX[n + 1];
         int const nYStart = VnLineGridY[n];
         int const nYEnd = VnLineGridY[n + 1];

         // Interpolate between pairs of points using a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm). Bresenham's algorithm gave occasional gaps
         double dXInc = nXEnd - nXStart;
         double dYInc = nYEnd - nYStart;
         double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

         dXInc /= dLength;
         dYInc /= dLength;

         double dX = nXStart;
         double dY = nYStart;

         // Process each interpolated point
         for (int m = 0; m <= nRound(dLength); m++)
         {
            int const nX = nRound(dX);
            int const nY = nRound(dY);

            // Have we hit a coastline cell?
            if (bIsWithinValidGrid(nX, nY) && m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               nCoastX = nX;
               nCoastY = nY;
               break;
            }

            // Two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
            if (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline())
            {
               nCoastX = nX;
               nCoastY = nY + 1;
               break;
            }

            // Increment for next time
            dX += dXInc;
            dY += dYInc;
         }

         // Did we find an intersection?
         if (nCoastX != -1)
            break;
      }

      // Did we find an intersection?
      if (nCoastX == -1)
      {
         // Nope
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": No line/coast intersection found" << endl;

         return RTN_ERR_SEDIMENT_INPUT_EVENT;
      }

      // OK we have an intersection of the line and coast. We will input the sediment here. Get landform and top layer
      CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform();
      int const nTopLayer = m_pRasterGrid->m_Cell[nCoastX][nCoastY].nGetNumOfTopLayerAboveBasement();

      // Is this intersection point in a polygon?
      int const nThisPoly = m_pRasterGrid->m_Cell[nCoastX][nCoastY].nGetPolygonID();
      int nThisPolyCoast = INT_NODATA;
      if (nThisPoly != INT_NODATA)
      {
         // Yes we are in a polygon, so get the coast ID of the polygon for this cell
         nThisPolyCoast = m_pRasterGrid->m_Cell[nCoastX][nCoastY].nGetPolygonCoastID();

         // Safety check
         if (nThisPolyCoast == INT_NODATA)
            return RTN_ERR_CELL_IN_POLY_BUT_NO_POLY_COAST;
      }

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         LogStream << m_ulIter << ": line/coast intersection is at [" << nCoastX << "][" << nCoastY << "] = {" << dGridXToExtCRSX(nCoastX) << ", " << dGridYToExtCRSY(nCoastY) << "}";

         if (nThisPoly != INT_NODATA)
            LogStream << " which is within coast " << nThisPolyCoast << " polygon " << nThisPoly << endl;
         else
            LogStream << " which is not within a polygon" << endl;
      }

      // Is some fine unconsolidated sediment being input?
      double const dFineDepth = dFineSedVol / m_dCellArea;
      if (dFineDepth > 0)
      {
         // Yes, so add to this cell's fine unconsolidated sediment
         m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFineSedimentInputDepth(dFineDepth);

         // And update the sediment top elevation value
         m_pRasterGrid->m_Cell[nCoastX][nCoastY].CalcAllLayerElevsAndD50();

         if (nThisPoly != INT_NODATA)
            // Add to this polygon's fine sediment input total
           m_VCoast[nThisPolyCoast].pGetPolygon(nThisPoly)->SetSedimentInputUnconsFine(dFineDepth);

         // Add to the this-iteration total of fine sediment input
         m_dThisiterUnconsFineInput += dFineDepth;

         // And assign the cell's landform category
         pLandform->SetLFCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
      }

      // Is some sand-sized unconsolidated sediment being input?
      double const dSandDepth = dSandSedVol / m_dCellArea;
      if (dSandDepth > 0)
      {
         // Yes, so add to this cell's sand unconsolidated sediment
         m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandSedimentInputDepth(dSandDepth);

         // And update the sediment top elevation value
         m_pRasterGrid->m_Cell[nCoastX][nCoastY].CalcAllLayerElevsAndD50();

         if (nThisPoly != INT_NODATA)
            // Add to this polygon's sand sediment input total
            m_VCoast[nThisPolyCoast].pGetPolygon(nThisPoly)->SetSedimentInputUnconsSand(dSandDepth);

         // Add to the this-iteration total of sand sediment input
         m_dThisiterUnconsSandInput += dSandDepth;

         // And assign the cell's landform category
         pLandform->SetLFCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
      }

      // Is some coarse unconsolidated sediment being input?
      double const dCoarseDepth = dCoarseSedVol / m_dCellArea;
      if (dCoarseDepth > 0)
      {
         // Yes, so add to this cell's coarse unconsolidated sediment
         m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseSedimentInputDepth(dCoarseDepth);

         // And update the sediment top elevation value
         m_pRasterGrid->m_Cell[nCoastX][nCoastY].CalcAllLayerElevsAndD50();

         if (nThisPoly != INT_NODATA)
            // Add to this polygon's coarse sediment input total
            m_VCoast[nThisPolyCoast].pGetPolygon(nThisPoly)->SetSedimentInputUnconsCoarse(dCoarseDepth);

         // Add to the this-iteration total of coarse sediment input
         m_dThisiterUnconsCoarseInput += dCoarseDepth;

         // And assign the cell's landform category
         pLandform->SetLFCategory(LF_SEDIMENT_INPUT_UNCONSOLIDATED);
      }

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << "; depth of fine sediment added = " << dFineDepth << " m, depth of sand sediment added = " << dSandDepth << " m, depth of coarse sediment added = " << dCoarseDepth << " m" << endl;
   }

   return RTN_OK;
}
