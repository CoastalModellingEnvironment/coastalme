/*!
   \file create_polygons.cpp
   \brief Creates coast polygons for sediment transport calcs
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
#include <assert.h>

#include <iostream>
using std::endl;

#include <stack>
using std::stack;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2d_point.h"
#include "2di_point.h"

//===============================================================================================================================
//! Create polygons, and mark the polygon boundaries on the raster grid
//===============================================================================================================================
int CSimulation::nCreateAllPolygons(void)
{
   LogStream << endl << m_ulIter << ": Creating polygons" << endl;

   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nNode = -1;
      int nNextProfile = -1;
      int nPolygon = -1;

      // Do this for every point on the coastline (except the last point)
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      for (int nCoastPoint = 0; nCoastPoint < nCoastSize - 1; nCoastPoint++)
      {
         if (! m_VCoast[nCoast].bIsProfileAtCoastPoint(nCoastPoint))
            continue;

         // OK, this coast point is the start of a coastline-normal profile
         CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nCoastPoint);

         if (pThisProfile->bOKIncStartAndEndOfCoast())
         {
            // This profile is OK, so we will start a polygon here and extend it down-coast (i.e. along the coast in the direction of increasing coastline point numbers)
            int const nThisProfile = pThisProfile->nGetProfileID();

            // This will be the coast ID number of the polygon, and also the polygon's along-coast sequence
            nPolygon++;

            // Now get a pointer to the next (down-coast) profile
            CGeomProfile* pNextProfile = pThisProfile->pGetDownCoastAdjacentProfile();

            bool bNextProfileIsOK = false;

            do
            {
               // if (pNextProfile == NULL)       // TODO THIS GIVES CPPCHECK ERROR
               // {
               // LogStream << m_ulIter << ": nThisProfile = " << nThisProfile << " invalid down-coast adjacent profile, trying next down-coast adjacent profile" << endl;
               //
               //    // Try the next-after-next profile
               // CGeomProfile* pNextNextProfile = pNextProfile->pGetDownCoastAdjacentProfile();
               // pNextProfile = pNextNextProfile;
               // continue;
               // }

               // Get the ID of the next (down-coast) profile
               nNextProfile = pNextProfile->nGetProfileID();

               // Is the next profile OK?
               bNextProfileIsOK = pNextProfile->bOKIncStartAndEndOfCoast();

               if (! bNextProfileIsOK)
               {
                  // Nope, the next profile is not OK
                  // LogStream << m_ulIter << ": down-coast adjacent profile = " << nNextProfile << " is not OK" << endl;

                  // So try the following profile
                  CGeomProfile* pNextNextProfile = pNextProfile->pGetDownCoastAdjacentProfile();
                  pNextProfile = pNextNextProfile;
               }

            } while (! bNextProfileIsOK);

            // LogStream << "Profile " << pNextProfile->nGetProfileID() << " is OK" << endl;

            // Get the coast point at which this next profile starts
            int const nNextProfileCoastPoint = pNextProfile->nGetCoastPoint();

            // Calculate half the along-coast distance (in coast points) between this profile and the next (i.e. down-coast) profile
            int const nDist = (nNextProfileCoastPoint - nCoastPoint) / 2;

            // OK, set the node point in the coast object. We do this now, instead of earlier on, since some profiles (i.e. polygon boundaries) may have been marked as invalid
            int const nNodePoint = nCoastPoint + nDist;
            m_VCoast[nCoast].SetPolygonNode(nNodePoint, ++nNode);

            // Get the grid CRS coordinates of the coast node
            CGeom2DIPoint const PtiNode = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nNodePoint);

            // Get some defaults (assuming for now that this polygon is not approximately triangular i.e. both normals do not meet)
            int nNextProfileEnd = pNextProfile->nGetProfileSize() - 1;
            int nThisProfileEnd = pThisProfile->nGetProfileSize() - 1;
            bool bMeetsAtAPoint = false;
            CGeom2DPoint PtCoastwardTip;

            //                // DEBUG CODE =============================================================================================
            // CGeom2DPoint PtNextEndTmp = *pNextProfile->pPtGetPointInProfile(nNextProfileEnd);
            // double dXTmp = dExtCRSXToGridX(PtNextEndTmp.dGetX());
            // double dYTmp = dExtCRSYToGridY(PtNextEndTmp.dGetY());
            // assert(bIsWithinValidGrid(dXTmp, dYTmp));
            //
            // CGeom2DPoint PtThisEndTmp = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
            // dXTmp = dExtCRSXToGridX(PtThisEndTmp.dGetX());
            // dYTmp = dExtCRSYToGridY(PtThisEndTmp.dGetY());
            // assert(bIsWithinValidGrid(dXTmp, dYTmp));
            //                // DEBUG CODE =============================================================================================

            // Now check to see if the two normals do meet i.e. if they are coincident
            if (pThisProfile->bFindProfileInCoincidentProfiles(nNextProfile))
            {
               // Yes they do meet
               bMeetsAtAPoint = true;
               int nTmpThisProfileEnd;
               int nTmpNextProfileEnd;

               // Find the most coastward point at which this normal and the previous normal touch. If they do not touch, the polygon requires a 'joining line'
               pThisProfile->GetMostCoastwardSharedLineSegment(nNextProfile, nTmpThisProfileEnd, nTmpNextProfileEnd);

               if (nTmpThisProfileEnd == -1)
               {
                  LogStream << m_ulIter << ": " << ERR << "profile " << nNextProfile << " should be coincident with profile " << nThisProfile << " but was not found" << endl;
                  return RTN_ERR_BAD_MULTILINE;
               }

               // Safety check: make sure that nThisProfileEnd is no bigger than pThisProfile->nGetProfileSize()-1, and the same for nNextProfileEnd
               nThisProfileEnd = tMin(nThisProfileEnd, nTmpThisProfileEnd);
               nNextProfileEnd = tMin(nNextProfileEnd, nTmpNextProfileEnd);

               PtCoastwardTip = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
            }

            // Create the vector in which to store the polygon's boundary (external CRS). Note that these points are not in sequence
            vector<CGeom2DPoint> PtVBoundary;

            // Start appending points: begin by appending the points in this normal, in reverse order
            for (int i = nThisProfileEnd; i >= 0; i--)
            {
               CGeom2DPoint const PtThis = *pThisProfile->pPtGetPointInProfile(i);
               PtVBoundary.push_back(PtThis);
            }

            // Next add coast points: from the start point of this normal, moving down-coast as far as the down-coast norma
            for (int i = nCoastPoint; i <= nNextProfileCoastPoint; i++)
            {
               CGeom2DPoint const PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(i);
               PtVBoundary.push_back(PtThis);
            }

            // Append the points in the down-coast (next) normal
            for (int i = 0; i <= nNextProfileEnd; i++)
            {
               CGeom2DPoint const PtThis = *pNextProfile->pPtGetPointInProfile(i);
               PtVBoundary.push_back(PtThis);
            }

            // Finally, append the end point of this normal
            CGeom2DPoint const PtThis = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
            PtVBoundary.push_back(PtThis);

            // Now identify the 'anti-node', this is the seaward point 'opposite' the polygon's coastal node
            CGeom2DIPoint PtiAntiNode;

            if (bMeetsAtAPoint)
               PtiAntiNode = PtiExtCRSToGridRound(&PtCoastwardTip);
            else
            {
               CGeom2DPoint const PtAvg = PtAverage(pThisProfile->pPtGetPointInProfile(nThisProfileEnd), pNextProfile->pPtGetPointInProfile(nNextProfileEnd));
               PtiAntiNode = PtiExtCRSToGridRound(&PtAvg);
            }

            // Is this a start-of-coast or an end-of-coast polygon?
            bool bStartCoast = false;
            bool bEndCoast = false;

            if (pThisProfile->bStartOfCoast())
               bStartCoast = true;

            if (pNextProfile->bEndOfCoast())
               bEndCoast = true;

            // Create the coast polygon object, store it, and get a pointer to it
            CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pPolyCreateAndAppendPolygon(nPolygon, nNodePoint, &PtiNode, &PtiAntiNode, nThisProfile, nNextProfile, &PtVBoundary, nThisProfileEnd + 1, nNextProfileEnd + 1, bStartCoast, bEndCoast);

            // Save the profile-end vertices (but omit the last one if the profiles meet at a point)
            pPolygon->AppendVertex(pThisProfile->pPtiGetStartPoint());
            pPolygon->AppendVertex(pThisProfile->pPtiGetEndPoint());
            pPolygon->AppendVertex(pNextProfile->pPtiGetStartPoint());

            if (! bMeetsAtAPoint)
               pPolygon->AppendVertex(pNextProfile->pPtiGetEndPoint());

            // // DEBUG CODE =================================================================================
            // LogStream << m_ulIter << ": vertices for coast " << nCoast << " polygon = " << nPolygon << endl;
            // for (int n = 0; n < pPolygon->nGetNumVertices(); n++)
            // LogStream << "[" << pPolygon->PtiGetVertex(n).nGetX() << "][" << pPolygon->PtiGetVertex(n).nGetY() << "]\t";
            // LogStream << endl;
            // // DEBUG CODE =================================================================================

            // assert(nPolygon < m_VCoast[nCoast].nGetNumPolygons());

            // Now rasterize the polygon boundaries: first, the coastline. This is necessary so that sand/coarse sediment derived from platform erosion of the coast cells is correctly added to the containing polygon's unconsolidated sediment
            for (int i = nCoastPoint; i <= nNextProfileCoastPoint; i++)
            {
               CGeom2DIPoint const PtiToMark = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(i);
               m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetCoastAndPolygonID(nCoast, nPolygon);
            }

            // Do the down-coast normal profile (does the whole length, including any shared line segments. So some cells are marked twice, however this is not a problem)
            int nCellsInProfile = pNextProfile->nGetNumCellsInProfile();

            for (int i = 0; i < nCellsInProfile; i++)
            {
               CGeom2DIPoint const PtiToMark = *pNextProfile->pPtiGetCellInProfile(i);
               m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetCoastAndPolygonID(nCoast, nPolygon);
            }

            // Do this normal profile (again does the whole length, including any shared line segments)
            nCellsInProfile = pThisProfile->nGetNumCellsInProfile();

            for (int i = 0; i < nCellsInProfile; i++)
            {
               CGeom2DIPoint const PtiToMark = *pThisProfile->pPtiGetCellInProfile(i);
               m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetCoastAndPolygonID(nCoast, nPolygon);
            }

            // If the polygon doesn't meet at a point at its seaward end, also need to rasterize the 'joining line'
            if (! bMeetsAtAPoint)
            {
               CGeom2DIPoint const PtiDownCoastNormalEnd = *pNextProfile->pPtiGetEndPoint(); // Grid CRS
               CGeom2DIPoint const PtiUpCoastNormalEnd = *pThisProfile->pPtiGetEndPoint();   // Grid CRS

               RasterizePolygonJoiningLine(nCoast, &PtiUpCoastNormalEnd, &PtiDownCoastNormalEnd, nPolygon);

               // pPolygon->SetNotPointed();
            }

            // }

            // nNextProfile = nThisProfile;
            nCoastPoint = nNextProfileCoastPoint - 1;

            if (bEndCoast)
               break;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Puts a polygon 'joining line' (the line which is the seaward boundary of the polygon, if the polygon doesn't meet at a point) onto the raster grid
//===============================================================================================================================
void CSimulation::RasterizePolygonJoiningLine(int nCoast, CGeom2DIPoint const* pPt1, CGeom2DIPoint const* pPt2, int const nPoly)
{
   // The start point of the line (grid CRS)
   int const nXStart = pPt1->nGetX();
   int const nYStart = pPt1->nGetY();

   // The end point of the line (grid CRS)
   int const nXEnd = pPt2->nGetX();
   int const nYEnd = pPt2->nGetY();

   // Safety check, in case the two points are identical (can happen due to rounding errors)
   if ((nXStart == nXEnd) && (nYStart == nYEnd))
      return;

   // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
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
      int nX = nRound(dX);
      int nY = nRound(dY);

      // Safety check
      if (! bIsWithinValidGrid(nX, nY))
         KeepWithinValidGrid(nXStart, nYStart, nX, nY);

      // assert(nPoly < m_VCoast[0].nGetNumPolygons());

      // Mark this point on the raster grid
      m_pRasterGrid->m_Cell[nX][nY].SetCoastAndPolygonID(nCoast, nPoly);

      // And increment for next time
      dX += dXInc;
      dY += dYInc;
   }
}

//===============================================================================================================================
//! Marks cells of the raster grid that are within each coastal polygon. The cell-by-cell fill (aka 'floodfill') code used here is adapted from an example by Lode Vandevenne (http://lodev.org/cgtutor/floodfill.html#Scanline_Floodfill_Algorithm_With_Stack)
//===============================================================================================================================
void CSimulation::MarkPolygonCells(void)
{
   // // DEBUG CODE =================================================================================
   // vector<CGeom2DIPoint> VPtiCentroids;
   // vector<int> VnID;
   //
   // // Do this for each coast
   // for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   // {
   //    // Do this for every coastal polygon
   // for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
   // {
   // CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
   //
   // int nPolyID = pPolygon->nGetPolygonCoastID();
   // VnID.push_back(nPolyID);
   //
   // CGeom2DIPoint PtiStart = pPolygon->PtiGetFillStartPoint();
   // VPtiCentroids.push_back(PtiStart);
   // }
   // }
   //
   // string strOutFile = m_strOutPath;
   // strOutFile += "00_polygon_flood_fill_start_point_";
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
   // pdRaster[nn] = INT_NODATA;
   //
   // for (int n = 0; n < static_cast<int>(VPtiCentroids.size()); n++)
   // {
   // if ((VPtiCentroids[n].nGetX() == nX) && (VPtiCentroids[n].nGetY() == nY))
   // {
   // pdRaster[nn] = VnID[n];
   // LogStream << m_ulIter << ": cell-by-cell fill start point for polygon " << VnID[n] << " is [" << nX << "][" << nY << "]" << endl;
   // break;
   // }
   // }
   // nn++;
   // }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   //
   // if (nRet == CE_Failure)
   // LogStream << nRet << endl;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Do this for every coastal polygon, in along-coast sequence
      int const nPolygons = m_VCoast[nCoast].nGetNumPolygons();

      for (int nPoly = 0; nPoly < nPolygons; nPoly++)
      {
         int nCellsInPolygon = 0;
         double dTotDepth = 0;
         double dStoredUnconsFine = 0;
         double dStoredUnconsSand = 0;
         double dStoredUnconsCoarse = 0;
         double dStoredConsFine = 0;
         double dStoredConsSand = 0;
         double dStoredConsCoarse = 0;
         double dSedimentInputFine = 0;
         double dSedimentInputSand = 0;
         double dSedimentInputCoarse = 0;

         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         int const nPolyID = pPolygon->nGetPolygonCoastID();

         // LogStream << m_ulIter << ": in MarkPolygonCells() nPoly = " << nPoly << " nPolyID = " << nPolyID << endl;

         // Create an empty stack
         stack<CGeom2DIPoint> PtiStack;

         //          // Since the polygon's vector boundary does not coincide exactly with the polygon's raster boundary, and the point-in-polygon check gives an indeterminate result if the point is exactly on the polygon's boundary, for safety we must construct a vector 'inner buffer' which is smaller than, and inside, the vector boundary TODO *** STILL NEEDED?
         // int nHand = m_VCoast[nCoast].nGetSeaHandedness();
         // int nSize = pPolygon->nGetBoundarySize();
         // vector<CGeom2DPoint> PtVInnerBuffer;
         //
         // for (int i = 0; i < nSize-1; i++)
         // {
         // int j = i + 1;
         // if (i == nSize-2)       // We must ignore the duplicated node point
         // j = 0;
         //
         // CGeom2DPoint PtThis = *pPolygon->pPtGetBoundaryPoint(i);
         // CGeom2DPoint PtNext = *pPolygon->pPtGetBoundaryPoint(j);
         //
         //             // Safety check
         // if (PtThis == PtNext)
         // continue;
         //
         // CGeom2DPoint PtBuffer = PtGetPerpendicular(&PtThis, &PtNext, m_dCellSide, nHand);
         //
         // PtVInnerBuffer.push_back(PtBuffer);
         // }

         // // DEBUG CODE ==============================================================================================
         // LogStream << endl << m_ulIter << ": coast " << nCoast << ", polygon " << nPoly << endl;
         // LogStream << "Boundary\t\t\tBuffer" << endl;
         // for (int i = 0; i < pPolygon->nGetBoundarySize()-1; i++)
         // LogStream << "{" << pPolygon->pPtGetBoundaryPoint(i)->dGetX() << ", " << pPolygon->pPtGetBoundaryPoint(i)->dGetY() << "}\t{" << PtVInnerBuffer[i].dGetX() << ", " << PtVInnerBuffer[i].dGetY() << "}" << endl;
         // LogStream << endl;
         // // DEBUG CODE ==============================================================================================

         // Use the centroid as the start point for the cell-by-cell fill procedure
         CGeom2DIPoint const PtiStart = pPolygon->PtiGetFillStartPoint();

         // // Is the centroid within the inner buffer?
         // if (! bIsWithinPolygon(&PtStart, &PtVInnerBuffer))
         // {
         //    // No, it is not: the polygon must be a concave polygon. So keep looking for a point which is definitely inside the polygon, using an alternative method
         // PtStart = PtFindPointInPolygon(&PtVInnerBuffer, pPolygon->nGetPointInPolygonSearchStartPoint());
         // }
         //
         // // Safety check (PtFindPointInPolygon() returns CGeom2DPoint(DBL_NODATA, DBL_NODATA) if it cannot find a valid start point)
         // if (bFPIsEqual(PtStart.dGetX(), DBL_NODATA, TOLERANCE))
         // {
         // LogStream << m_ulIter << ": " << ERR << "could not find a cell-by-cell fill start point for coast " << nCoast << ", polygon " << nPoly << endl;
         // break;
         // }

         // We have a cell-by-cell fill start point which is definitely within the polygon so push this point onto the stack
         // CGeom2DIPoint PtiStart;
         // PtiStart.SetX(nRound(PtStart.dGetX()));               // Grid CRS
         // PtiStart.SetY(nRound(PtStart.dGetY()));               // Grid CRS

         PtiStack.push(PtiStart);

         // LogStream << m_ulIter << ": filling polygon " << nPoly << " from [" << PtiStart.nGetX() << "][" << PtiStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiStart.nGetY()) << "}" << endl;

         // Then do the cell-by-cell fill: loop until there are no more cell coordinates on the stack
         while (! PtiStack.empty())
         {
            CGeom2DIPoint const Pti = PtiStack.top();
            PtiStack.pop();

            int nX = Pti.nGetX();
            int const nY = Pti.nGetY();

            while ((nX >= 0) && (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() == INT_NODATA))
               nX--;

            nX++;

            bool bSpanAbove = false;
            bool bSpanBelow = false;

            while ((nX < m_nXGridSize) && (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() == INT_NODATA))
            {
               // assert(nPolyID < m_VCoast[nCoast].nGetNumPolygons());

               // Mark the cell as being in this polygon and this coast
               m_pRasterGrid->m_Cell[nX][nY].SetCoastAndPolygonID(nCoast, nPolyID);

               // LogStream << "Marked [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

               // Increment the running totals for this polygon
               // dCliffCollapseErosionFine += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionFine();
               // dCliffCollapseErosionSand += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionSand();
               // dCliffCollapseErosionCoarse += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionCoarse();
               // dCliffCollapseTalusSand += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseSandTalusDeposition();
               // dCliffCollapseTalusCoarse += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseCoarseTalusDeposition();

               // Get the number of the highest layer with non-zero thickness
               int const nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

               // Safety check
               if ((nThisLayer != NO_NONZERO_THICKNESS_LAYERS) && (nThisLayer != INT_NODATA))
               {
                  // Increment some more running totals for this polygon TODO 066 should this be for ALL layers above the basement?
                  dStoredUnconsFine += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
                  dStoredUnconsSand += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
                  dStoredUnconsCoarse += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();

                  dStoredConsFine += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetFineDepth();
                  dStoredConsSand += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetSandDepth();
                  dStoredConsCoarse += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();

                  // Add to the start-iteration total of suspended fine sediment within polygons
                  m_dStartIterSuspFineInPolygons += m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();

                  // Add to the total of sediment derived from sediment input events
                  dSedimentInputFine += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFineSedimentInputDepth();
                  dSedimentInputSand += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSandSedimentInputDepth();
                  dSedimentInputCoarse += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarseSedimentInputDepth();
               }

               nCellsInPolygon++;
               dTotDepth += m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();

               if ((! bSpanAbove) && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY - 1].nGetPolygonID() == INT_NODATA))
               {
                  PtiStack.push(CGeom2DIPoint(nX, nY - 1));
                  bSpanAbove = true;
               }

               else if (bSpanAbove && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY - 1].nGetPolygonID() != INT_NODATA))
               {
                  bSpanAbove = false;
               }

               if ((! bSpanBelow) && (nY < m_nYGridSize - 1) && (m_pRasterGrid->m_Cell[nX][nY + 1].nGetPolygonID() == INT_NODATA))
               {
                  PtiStack.push(CGeom2DIPoint(nX, nY + 1));
                  bSpanBelow = true;
               }

               else if (bSpanBelow && (nY < m_nYGridSize - 1) && (m_pRasterGrid->m_Cell[nX][nY + 1].nGetPolygonID() != INT_NODATA))
               {
                  bSpanBelow = false;
               }

               nX++;
            }
         }

         // Store this polygon's stored unconsolidated sediment depths
         pPolygon->SetPreExistingUnconsFine(dStoredUnconsFine);
         pPolygon->SetPreExistingUnconsSand(dStoredUnconsSand);
         pPolygon->SetPreExistingUnconsCoarse(dStoredUnconsCoarse);

         // Store this polygon's values for unconsolidated sediment derived from sediment input event(s)
         pPolygon->SetSedimentInputUnconsFine(dSedimentInputFine);
         pPolygon->SetSedimentInputUnconsSand(dSedimentInputSand);
         pPolygon->SetSedimentInputUnconsCoarse(dSedimentInputCoarse);

         // Store this polygon's stored consolidated sediment depths
         pPolygon->SetPreExistingConsFine(dStoredConsFine);
         pPolygon->SetPreExistingConsSand(dStoredConsSand);
         pPolygon->SetPreExistingConsCoarse(dStoredConsCoarse);

         // Store the number of cells in the interior of the polygon
         pPolygon->SetNumCellsInPolygon(nCellsInPolygon);
         // LogStream << m_ulIter << ": in MarkPolygonCells() N cells = " << nCellsInPolygon << " in polygon " << nPolyID << endl;

         // Calculate the total volume of seawater on the polygon (m3) and store it
         double const dSeaVolume = dTotDepth * m_dCellSide;
         pPolygon->SetSeawaterVolume(dSeaVolume);
      }
   }

   // // DEBUG CODE ===========================================================================================================
   // if (m_ulIter == 109)
   // {
   // string strOutFile = m_strOutPath + "00_polygon_fill_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   //
   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   //
   // int nNumPoly = m_VCoast[0].nGetNumPolygons();
   // vector<int> nVPerPoly(nNumPoly, 0);
   // int nNotInPoly = 0;
   //
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   // for (int nX = 0; nX < m_nXGridSize; nX++)
   // {
   // int nID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
   // if (nID == INT_NODATA)
   // nNotInPoly++;
   // else
   // nVPerPoly[nID]++;
   //
   // pdRaster[n++] = nID;
   // }
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
   //
   // for (int nn = 0; nn < nNumPoly; nn++)
   // LogStream << m_ulIter << ": polygon " << nn << " has " << nVPerPoly[nn] << " cells" << endl;
   // LogStream << m_ulIter << " Number of cells not in any polygon = " << nNotInPoly << endl;
   // LogStream << "==================" << endl;
   // }
   // // DEBUG CODE ===========================================================================================================
}

//===============================================================================================================================
//! For between-polygon potential sediment routing: find which are the adjacent polygons, and calc the length of the shared normal between this polygon and the adjacent polygons TODO 012 Will need to change this when length of coastline-normal profiles (and so polygon seaward length) is determined by depth of closure
//===============================================================================================================================
int CSimulation::nDoPolygonSharedBoundaries(void)
{
   // Do this for every coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

      // // DEBUG CODE =================
      // for (int m = 0; m < m_VCoast[nCoast].nGetNumProfiles(); m++)
      // {
      // CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(m);
      //
      // LogStream << m << "\t" << pProfile->nGetProfileID() << "\t";
      //
      // int nPointsInProfile = pProfile->nGetProfileSize();
      //
      // for (int nPoint = 0; nPoint < nPointsInProfile; nPoint++)
      // {
      // CGeom2DPoint Pt = *pProfile->pPtGetPointInProfile(nPoint);
      // LogStream << "{" << Pt.dGetX() << ", " << Pt.dGetY() << "}";
      // }
      // LogStream << endl;
      // }
      // // DEBUG CODE =================

      // Do this for every coastal polygon, in along-coast sequence
      for (int nn = 0; nn < nNumPolygons; nn++)
      {
         CGeomCoastPolygon* pThisPolygon = m_VCoast[nCoast].pGetPolygon(nn);
         int const nThisPolygon = pThisPolygon->nGetPolygonCoastID();

         vector<int> nVUpCoastAdjacentPolygon;
         vector<int> nVDownCoastAdjacentPolygon;

         vector<double> dVUpCoastBoundaryShare;
         vector<double> dVDownCoastBoundaryShare;

         // First deal with down-coast adjacent polygons
         if (pThisPolygon->bIsCoastEndPolygon())
         {
            // We are at the end of the coastline, no other polygon is adjacent to the down-coast profile of the end-of-coast polygon
            nVDownCoastAdjacentPolygon.push_back(INT_NODATA);
            dVDownCoastBoundaryShare.push_back(1);

            // Store in the polygon
            pThisPolygon->SetDownCoastAdjacentPolygons(&nVDownCoastAdjacentPolygon);
            pThisPolygon->SetDownCoastAdjacentPolygonBoundaryShares(&dVDownCoastBoundaryShare);
         }
         else
         {
            // We are not at the end of the coastline, so there is at least one other polygon adjacent to the down-coast profile of this polygon
            int const nDownCoastProfile = pThisPolygon->nGetDownCoastProfile();
            int const nPointsInProfile = pThisPolygon->nGetNumPointsUsedDownCoastProfile();

            double dDownCoastTotBoundaryLen = 0;

            CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nDownCoastProfile);

            for (int nPoint = 0; nPoint < nPointsInProfile - 1; nPoint++)
            {
               CGeom2DPoint const PtStart = *pProfile->pPtGetPointInProfile(nPoint);
               CGeom2DPoint const PtEnd = *pProfile->pPtGetPointInProfile(nPoint + 1);

               // Calculate the length of this segment of the normal profile. Note that it should not be zero, since we checked for duplicate points when creating profiles
               double const dDistBetween = dGetDistanceBetween(&PtStart, &PtEnd);

               // Find out which polygon is adjacent to each line segment of the polygon's down-coast profile boundary. The basic approach used is to count the number of coincident profiles in each line segment, and (because we are going down-coast) add this number to 'this' polygon's number. However, some of these coincident profiles may be invalid, so we must count only the valid co-incident profiles
               int const nNumCoinc = pProfile->nGetNumCoincidentProfilesInLineSegment(nPoint);

               // Safety check
               if (nNumCoinc < 0)
                  continue;

               int nNumValidCoinc = 0;

               for (int nCoinc = 0; nCoinc < nNumCoinc; nCoinc++)
               {
                  int const nProf = pProfile->nGetCoincidentProfileForLineSegment(nPoint, nCoinc);

                  // Safety check
                  if (nProf == -1)
                     continue;

                  CGeomProfile const* pProf = m_VCoast[nCoast].pGetProfile(nProf);

                  if (pProf->bProfileOK())
                     nNumValidCoinc++;
               }

               // First stab at calculating the number of the adjacent polygon
               int nAdj = nThisPolygon + nNumValidCoinc;

               // However, if 'this' polygon is close to the end of the coastline, we get polygon numbers greater than the number of polygons i.e. beyond the end of the coastline. If this happens, set the adjacent polygon to 'off-edge'
               if (nAdj >= nNumPolygons)
                  nAdj = INT_NODATA;

               // Safety check: is this adjacent polygon already present in the list of down-coast adjacent polygons?
               if (pThisPolygon->bDownCoastIsAlreadyPresent(nAdj))
                  continue;

               // Not already present
               nVDownCoastAdjacentPolygon.push_back(nAdj);

               dDownCoastTotBoundaryLen += dDistBetween;
               dVDownCoastBoundaryShare.push_back(dDistBetween);
            }

            // Calculate the down-coast boundary share
            for (unsigned int n = 0; n < dVDownCoastBoundaryShare.size(); n++)
            {
               // Safety check
               if (bFPIsEqual(dDownCoastTotBoundaryLen, 0.0, TOLERANCE))
               {
                  nVDownCoastAdjacentPolygon.pop_back();
                  dVDownCoastBoundaryShare.pop_back();
                  continue;
               }

               dVDownCoastBoundaryShare[n] /= dDownCoastTotBoundaryLen;
            }

            // Store in the polygon
            pThisPolygon->SetDownCoastAdjacentPolygons(&nVDownCoastAdjacentPolygon);
            pThisPolygon->SetDownCoastAdjacentPolygonBoundaryShares(&dVDownCoastBoundaryShare);
         }

         // Now deal with up-coast adjacent polygons
         if (pThisPolygon->bIsCoastStartPolygon())
         {
            // We are at the start of the coastline, no other polygon is adjacent to the up-coast profile of the start-of-coast polygon
            nVUpCoastAdjacentPolygon.push_back(INT_NODATA);
            dVUpCoastBoundaryShare.push_back(1);

            // Store in the polygon
            pThisPolygon->SetUpCoastAdjacentPolygons(&nVUpCoastAdjacentPolygon);
            pThisPolygon->SetUpCoastAdjacentPolygonBoundaryShares(&dVUpCoastBoundaryShare);
         }

         else
         {
            // We are not at the start of the coastline, so there is at least one other polygon adjacent to the up-coast profile of this polygon
            int const nUpCoastProfile = pThisPolygon->nGetUpCoastProfile();
            int const nPointsInProfile = pThisPolygon->nGetNumPointsUsedUpCoastProfile();

            double dUpCoastTotBoundaryLen = 0;

            CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nUpCoastProfile);

            for (int nPoint = 0; nPoint < nPointsInProfile - 1; nPoint++)
            {
               CGeom2DPoint const PtStart = *pProfile->pPtGetPointInProfile(nPoint);
               CGeom2DPoint const PtEnd = *pProfile->pPtGetPointInProfile(nPoint + 1);

               // Safety check
               if (PtStart == PtEnd)
                  // Should never get here, since we checked for duplicate points when creating profiles
                  continue;

               // Calculate the length of this segment of the normal profile
               double const dDistBetween = dGetDistanceBetween(&PtStart, &PtEnd);

               // Find out which polygon is adjacent to each line segment of the polygon's up-coast profile boundary. The basic approach used is to count the number of coincident profiles in each line segment, and (because we are going up-coast) subtract this number from 'this' polygon's number. However, some of these coincident profiles may be invalid, so we must count only the valid co-incident profiles
               int const nNumCoinc = pProfile->nGetNumCoincidentProfilesInLineSegment(nPoint);

               // Safety check
               if (nNumCoinc < 0)
                  continue;

               int nNumValidCoinc = 0;

               for (int nCoinc = 0; nCoinc < nNumCoinc; nCoinc++)
               {
                  int const nProf = pProfile->nGetCoincidentProfileForLineSegment(nPoint, nCoinc);

                  // Safety check
                  if (nProf == -1)
                     continue;

                  CGeomProfile const* pProf = m_VCoast[nCoast].pGetProfile(nProf);

                  if (pProf->bProfileOK())
                     nNumValidCoinc++;
               }

               // First stab at calculating the number of the adjacent polygon
               int nAdj = nThisPolygon - nNumValidCoinc;

               // However, if 'this' polygon is close to the start of the coastline, we get polygon numbers below zero i.e. beyond the start of the coastline. If this happens, set the adjacent polygon to 'off-edge'
               if (nAdj < 0)
                  nAdj = INT_NODATA;

               // Safety check: is this adjacent polygon already present in the list of up-coast adjacent polygons?
               if (pThisPolygon->bUpCoastIsAlreadyPresent(nAdj))
                  continue;

               // Not already present
               nVUpCoastAdjacentPolygon.push_back(nAdj);

               dUpCoastTotBoundaryLen += dDistBetween;
               dVUpCoastBoundaryShare.push_back(dDistBetween);
            }

            // Calculate the up-coast boundary share
            for (unsigned int n = 0; n < dVUpCoastBoundaryShare.size(); n++)
            {
               // Safety check
               if (bFPIsEqual(dUpCoastTotBoundaryLen, 0.0, TOLERANCE))
               {
                  nVUpCoastAdjacentPolygon.pop_back();
                  dVUpCoastBoundaryShare.pop_back();
                  continue;
               }

               dVUpCoastBoundaryShare[n] /= dUpCoastTotBoundaryLen;
            }

            // Store in the polygon
            pThisPolygon->SetUpCoastAdjacentPolygons(&nVUpCoastAdjacentPolygon);
            pThisPolygon->SetUpCoastAdjacentPolygonBoundaryShares(&dVUpCoastBoundaryShare);
         }

         // Finally, calculate the distance between the coast node and the antinode of the polygon
         double const dPolygonSeawardLen = dGetDistanceBetween(pThisPolygon->pPtiGetNode(), pThisPolygon->pPtiGetAntiNode());

         // And store it
         pThisPolygon->SetLength(dPolygonSeawardLen);

         // // DEBUG CODE ======================================================================================================================
         // assert(dVUpCoastBoundaryShare.size() == nVUpCoastAdjacentPolygon.size());
         // assert(dVDownCoastBoundaryShare.size() == nVDownCoastAdjacentPolygon.size());
         //
         // LogStream << m_ulIter << ": polygon = " << nPoly << (pPolygon->bIsPointed() ? " IS TRIANGULAR" : "") << endl;
         // LogStream << m_ulIter << ": coast " << nCoast << " polygon " << nThisPolygon << endl;
         //
         // int nUpCoastProfile = pThisPolygon->nGetUpCoastProfile();
         // CGeomProfile* pUpCoastProfile = m_VCoast[0].pGetProfile(nUpCoastProfile);
         // int nUpCoastProfileCells = pUpCoastProfile->nGetNumCellsInProfile();
         //
         // int nDownCoastProfile = pThisPolygon->nGetDownCoastProfile();
         // CGeomProfile* pDownCoastProfile = m_VCoast[0].pGetProfile(nDownCoastProfile);
         // int nDownCoastProfileCells = pDownCoastProfile->nGetNumCellsInProfile();
         //
         // LogStream << "\tUp-coast profile = " << nUpCoastProfile << " down-coast profile = " << nDownCoastProfile << endl;
         // LogStream << "\tN cells in up-coast profile = " << nUpCoastProfileCells << " N cells in down-coast profile = " << nDownCoastProfileCells << endl;
         //
         // LogStream << "\tThere are " << nVUpCoastAdjacentPolygon.size() << " UP-COAST adjacent polygon(s) = ";
         // for (unsigned int n = 0; n < nVUpCoastAdjacentPolygon.size(); n++)
         // LogStream << nVUpCoastAdjacentPolygon[n] << " ";
         // LogStream << endl;
         //
         // LogStream << "\tThere are " << nVDownCoastAdjacentPolygon.size() << " DOWN-COAST adjacent polygon(s) = ";
         // for (unsigned int n = 0; n < nVDownCoastAdjacentPolygon.size(); n++)
         // LogStream << nVDownCoastAdjacentPolygon[n] << " ";
         // LogStream << endl;
         //
         // double dUpCoastTotBoundaryLen = 0;
         // LogStream << "\tUP-COAST boundary share(s) = ";
         // for (unsigned int n = 0; n < dVUpCoastBoundaryShare.size(); n++)
         // {
         // dUpCoastTotBoundaryLen += dVUpCoastBoundaryShare[n];
         // LogStream << dVUpCoastBoundaryShare[n] << " ";
         // }
         // LogStream << endl;
         // LogStream << "\tTotal UP-COAST boundary length = " << dUpCoastTotBoundaryLen << endl;
         //
         // double dDownCoastTotBoundaryLen = 0;
         // LogStream << "\tDOWN-COAST boundary share(s) = ";
         // for (unsigned int n = 0; n < dVDownCoastBoundaryShare.size(); n++)
         // {
         // dDownCoastTotBoundaryLen += dVDownCoastBoundaryShare[n];
         // LogStream << dVDownCoastBoundaryShare[n] << " ";
         // }
         // LogStream << endl;
         // LogStream << "\tTotal DOWN-COAST boundary length = " << dDownCoastTotBoundaryLen << endl;
         // // DEBUG CODE ======================================================================================================================
      }
   }

   // // DEBUG CODE =================================================================================
   // for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   // {
   // int nNumPoly = m_VCoast[nCoast].nGetNumPolygons();
   // for (int nPoly = 0; nPoly < nNumPoly; nPoly++)
   // {
   // CGeomCoastPolygon const* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
   //
   // int nNumUpCoastPolygons = pPolygon->nGetNumUpCoastAdjacentPolygons();
   // LogStream << "Polygon " << nPoly << " up-coast polygon(s) = ";
   // for (int nAdj = 0; nAdj < nNumUpCoastPolygons; nAdj++)
   // {
   // int nAdjPoly = pPolygon->nGetUpCoastAdjacentPolygon(nAdj);
   // LogStream << nAdjPoly;
   // if (nAdjPoly > nNumPoly-1)
   // LogStream << "***";
   // }
   //
   // int nNumDownCoastPolygons = pPolygon->nGetNumDownCoastAdjacentPolygons();
   // LogStream << " down-coast polygon(s) = ";
   // for (int nAdj = 0; nAdj < nNumDownCoastPolygons; nAdj++)
   // {
   // int nAdjPoly = pPolygon->nGetDownCoastAdjacentPolygon(nAdj);
   // LogStream << nAdjPoly;
   // if (nAdjPoly > nNumPoly-1)
   // LogStream << "***";
   // }
   // LogStream << endl;
   // }
   // }
   // // DEBUG CODE =================================================================================

   return RTN_OK;
}

// //===============================================================================================================================
// //! Determines whether a point is within a polygon: however if the point is exactly on the edge of the polygon, then the result is indeterminate. Modified from code at http://alienryderflex.com/polygon/, our thanks to Darel Rex Finley (DarelRex@gmail.com)
// //===============================================================================================================================
// bool CSimulation::bIsWithinPolygon(CGeom2DPoint const* pPtStart, vector<CGeom2DPoint> const* pPtPoints)
// {
//    bool bOddNodes = false;
//
//    int const nPolyCorners = static_cast<int>(pPtPoints->size());
//    int j = nPolyCorners - 1;
//
//    double const dX = pPtStart->dGetX();
//    double const dY = pPtStart->dGetY();
//
//    for (int i = 0; i < nPolyCorners; i++)
//    {
//       double const dCorneriX = pPtPoints->at(i).dGetX();
//       double const dCorneriY = pPtPoints->at(i).dGetY();
//       double const dCornerjX = pPtPoints->at(j).dGetX();
//       double const dCornerjY = pPtPoints->at(j).dGetY();
//
//       if ((dCorneriY < dY && dCornerjY >= dY) || (dCornerjY < dY && dCorneriY >= dY))
//       {
//          if (dCorneriX + (dY - dCorneriY) / (dCornerjY - dCorneriY) * (dCornerjX - dCorneriX) < dX)
//          {
//             bOddNodes = ! bOddNodes;
//          }
//       }
//
//       j = i;
//    }
//
//    return bOddNodes;
// }

// //===============================================================================================================================
// //! Finds a point in a polygon: is guaranteed to succeed, as every strictly closed polygon has at least one triangle that is completely contained within the polygon. Derived from an algorithm at http://stackoverflow.com/questions/9797448/get-a-point-inside-the-polygon
// //===============================================================================================================================
// CGeom2DPoint CSimulation::PtFindPointInPolygon(vector<CGeom2DPoint> const* pPtPoints, int const nStartPoint)
// {
//    int const nPolySize = static_cast<int>(pPtPoints->size());
//    int nOffSet = 0;
//    CGeom2DPoint PtStart;
//
//    do
//    {
//       // Choose three consecutive points from the polygon
//       vector<CGeom2DPoint> nVTestPoints;
//
//       for (int n = 0; n < 3; n++)
//       {
//          int nIndex = n + nStartPoint + nOffSet;
//
//          if (nIndex > nPolySize - 1)
//             nIndex -= nPolySize;
//
//          // Safety check
//          if (nIndex < 0)
//             return CGeom2DPoint(DBL_NODATA, DBL_NODATA);
//
//          nVTestPoints.push_back(pPtPoints->at(nIndex));
//       }
//
//       // Increment ready for next time
//       nOffSet++;
//
//       // Safety check
//       if (nOffSet >= (nPolySize + 3))
//          return CGeom2DPoint(DBL_NODATA, DBL_NODATA);
//
//       // Check if the halfway point between the first and the third point is inside the polygon
//       PtStart = PtAverage(&nVTestPoints[0], &nVTestPoints[2]);
//    } while (! bIsWithinPolygon(&PtStart, pPtPoints));
//
//    return PtStart;
// }
