/*!
   \file gis_utils.cpp
   \brief Various GIS-related functions, requires GDAL
   \details Note re. coordinate systems used

   1. In the raster CRS, cell[0][0] is at the top left (NW) corner of the grid. Raster grid co-oordinate [0][0] is actually the top left (NW) corner of this cell.

   2. We assume that the grid CRS and external CRS have parallel axes. If they have not, see http://www.gdal.org/classGDALDataset.html which says that:

   To convert between pixel/line (P,L) raster space, and projection coordinates (Xp,Yp) space Xp = padfTransform[0] + padfTransform[1] + padfTransform[2]; Yp
   = padfTransform[3] + padfTransform[4] + padfTransform[5];

   In a north-up image, padfTransform[1] is the pixel width, and padfTransform[5] is the pixel height. The upper left corner of the upper left pixel is at position (padfTransform[0], padfTransform[3]).

   3. Usually, raster grid CRS values are integer, i.e. they refer to a point which is at the centroid of a cell. They may also be -ve or greater than m_nXGridSize-1 i.e. may refer to a point which lies outside any cell of the raster grid.

   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment. CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <assert.h>

#include <cstdio>

#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include <cmath>
using std::atan2;

#include <cfloat>
#include <climits>
#include <cstdint>

#include <cstring>
using std::strstr;

#include <gdal.h>
#include <gdal_priv.h>
#include <cpl_string.h>

#include "cme.h"
#include "coast.h"
#include "raster_grid.h"
#include "2d_point.h"
#include "2di_point.h"

//===============================================================================================================================
//! Given the integer X-axis ordinate of a cell in the raster grid CRS, returns the external CRS X-axis ordinate of the cell's centroid
//===============================================================================================================================
double CSimulation::dGridCentroidXToExtCRSX(int const nGridX) const
{
   // TODO 064
   return (m_dGeoTransform[0] + (nGridX * m_dGeoTransform[1]) + (m_dGeoTransform[1] / 2));
}

//===============================================================================================================================
//! Given the integer Y-axis ordinate of a cell in the raster grid CRS, returns the external CRS Y-axis ordinate of the cell's centroid
//===============================================================================================================================
double CSimulation::dGridCentroidYToExtCRSY(int const nGridY) const
{
   // TODO 064
   return (m_dGeoTransform[3] + (nGridY * m_dGeoTransform[5]) + (m_dGeoTransform[5] / 2));
}

//===============================================================================================================================
//! Transforms a pointer to a CGeom2DIPoint in the raster grid CRS (assumed to be the centroid of a cell) to the equivalent CGeom2DPoint in the external CRS
//===============================================================================================================================
CGeom2DPoint CSimulation::PtGridCentroidToExt(CGeom2DIPoint const *pPtiIn) const
{
   // TODO 064
   double const dX = m_dGeoTransform[0] + (pPtiIn->nGetX() * m_dGeoTransform[1]) + (m_dGeoTransform[1] / 2);
   double const dY = m_dGeoTransform[3] + (pPtiIn->nGetY() * m_dGeoTransform[5]) + (m_dGeoTransform[5] / 2);

   return CGeom2DPoint(dX, dY);
}

//===============================================================================================================================
//! Given a real-valued X-axis ordinate in the raster grid CRS (i.e. not the centroid of a cell), returns the external CRS X-axis ordinate
//===============================================================================================================================
double CSimulation::dGridXToExtCRSX(double const dGridX) const
{
   // TODO 064 Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
   return m_dGeoTransform[0] + (dGridX * m_dGeoTransform[1]) - 1;
}

//===============================================================================================================================
//! Given a real-valued Y-axis ordinate in the raster grid CRS (i.e. not the centroid of a cell), returns the external CRS Y-axis ordinate
//===============================================================================================================================
double CSimulation::dGridYToExtCRSY(double const dGridY) const
{
   // TODO 064 Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
   return m_dGeoTransform[3] + (dGridY * m_dGeoTransform[5]) - 1;
}

//===============================================================================================================================
//! Transforms an X-axis ordinate in the external CRS to the equivalent X-axis ordinate in the raster grid CRS (the result is not rounded, and so may not be integer, and may be outside the grid)
//===============================================================================================================================
double CSimulation::dExtCRSXToGridX(double const dExtCRSX) const
{
   // TODO 064
   return ((dExtCRSX - m_dGeoTransform[0]) / m_dGeoTransform[1]) - 1;
}

//===============================================================================================================================
//! Transforms a Y-axis ordinate in the external CRS to the equivalent Y-axis ordinate in the raster grid CRS (the result is not rounded, and so may not be integer, and may be outside the grid)
//===============================================================================================================================
double CSimulation::dExtCRSYToGridY(double const dExtCRSY) const
{
   // TODO 064
   return ((dExtCRSY - m_dGeoTransform[3]) / m_dGeoTransform[5]) - 1;
}

//===============================================================================================================================
//! Transforms a pointer to a CGeom2DPoint in the external CRS to the equivalent CGeom2DIPoint in the raster grid CRS (both values rounded). The result may be outside the grid
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiExtCRSToGridRound(CGeom2DPoint const* pPtIn) const
{
   // TODO 064
   int const nX = nRound(((pPtIn->dGetX() - m_dGeoTransform[0]) / m_dGeoTransform[1]) - 1);
   int const nY = nRound(((pPtIn->dGetY() - m_dGeoTransform[3]) / m_dGeoTransform[5]) - 1);

   return CGeom2DIPoint(nX, nY);
}

//===============================================================================================================================
//! Returns the distance (in external CRS) between two points
//===============================================================================================================================
double CSimulation::dGetDistanceBetween(CGeom2DPoint const* Pt1, CGeom2DPoint const* Pt2)
{
   double const dXDist = Pt1->dGetX() - Pt2->dGetX();
   double const dYDist = Pt1->dGetY() - Pt2->dGetY();

   return hypot(dXDist, dYDist);
}

//===============================================================================================================================
//! Returns the distance (in external CRS) between two points
//===============================================================================================================================
double CSimulation::dGetDistanceBetween(CGeom2DIPoint const* Pti1, CGeom2DIPoint const* Pti2)
{
   double const dXDist = Pti1->nGetX() - Pti2->nGetX();
   double const dYDist = Pti1->nGetY() - Pti2->nGetY();

   return hypot(dXDist, dYDist);
}

//===============================================================================================================================
//! Returns the distance (in external CRS) between two points
//===============================================================================================================================
double CSimulation::dGetDistanceBetween(double const dX1, double const dY1, double const dX2, double const dY2)
{
   double const dXDist = dX1 - dX2;
   double const dYDist = dY1 - dY2;

   return hypot(dXDist, dYDist);
}

//===============================================================================================================================
//! Returns twice the signed area of a triangle, defined by three points
//===============================================================================================================================
double CSimulation::dTriangleAreax2(CGeom2DPoint const* pPtA, CGeom2DPoint const* pPtB, CGeom2DPoint const* pPtC)
{
   return (pPtB->dGetX() - pPtA->dGetX()) * (pPtC->dGetY() - pPtA->dGetY()) - (pPtB->dGetY() - pPtA->dGetY()) * (pPtC->dGetX() - pPtA->dGetX());
}

//===============================================================================================================================
//! Checks whether the supplied point (an x-y pair, in the grid CRS) is within the raster grid, and is a valid cell (i.e. the basement DEM is not NODATA)
//===============================================================================================================================
bool CSimulation::bIsWithinValidGrid(int const nX, int const nY) const
{
   if ((nX < 0) || (nX >= m_nXGridSize))
      return false;

   if ((nY < 0) || (nY >= m_nYGridSize))
      return false;

   if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
      return false;

   return true;
}

//===============================================================================================================================
//! Checks whether the supplied point (a reference to a CGeom2DIPoint, in the grid CRS) is within the raster grid, and is a valid cell (i.e. the basement DEM is not NODATA)
//===============================================================================================================================
bool CSimulation::bIsWithinValidGrid(CGeom2DIPoint const* Pti) const
{
   int const nX = Pti->nGetX();
   int const nY = Pti->nGetY();

   return bIsWithinValidGrid(nX, nY);
}

//===============================================================================================================================
//! Constrains the supplied point (in the grid CRS) to be a valid cell within the raster grid
//===============================================================================================================================
void CSimulation::KeepWithinValidGrid(int& nX, int& nY) const
{
   nX = tMax(nX, 0);
   nX = tMin(nX, m_nXGridSize - 1);

   nY = tMax(nY, 0);
   nY = tMin(nY, m_nYGridSize - 1);
}

//===============================================================================================================================
//! Constrains the second supplied point (both are CGeom2DIPoints, in the grid CRS) to be a valid cell within the raster grid
//===============================================================================================================================
void CSimulation::KeepWithinValidGrid(CGeom2DIPoint const* Pti0, CGeom2DIPoint* Pti1) const
{
   KeepWithinValidGrid(Pti0->nGetX(), Pti0->nGetY(), *Pti1->pnGetX(), *Pti1->pnGetY());
}

//===============================================================================================================================
//! Given the X and Y co-ordinates of a point in the ext CRS, calculates the grid CRS equivalent: both are constrained (in necessary) to be within the valid grid
//===============================================================================================================================
void CSimulation::KeepWithinValidGrid(double& dX, double& dY, int& nX, int& nY) const
{
   nX = nRound(dExtCRSXToGridX(dX));
   nY = nRound(dExtCRSYToGridY(dY));

   if (bIsWithinValidGrid(nX, nY))
      return;

   KeepWithinValidGrid(nX, nY);

   dX = dGridCentroidXToExtCRSX(nX);
   dY = dGridCentroidYToExtCRSY(nY);
}

//===============================================================================================================================
//! Given two points in the grid CRS (the points assumed not to be coincident), this routine modifies the value of the second point so that it is on a line joining the original two points and is a valid cell within the raster grid. However in some cases (e.g. if the first point is at the edge of the valid part of the raster grid) then the second cell will be coincident with the
//! first cell, and the line joining them is thus of zero length. The calling routine has to be able to handle this
//===============================================================================================================================
void CSimulation::KeepWithinValidGrid(int nX0, int nY0, int& nX1, int& nY1) const
{
   // Safety check: make sure that the first point is within the valid grid
   if (nX0 >= m_nXGridSize)
      nX0 = m_nXGridSize - 1;

   else if (nX0 < 0)
      nX0 = 0;

   if (nY0 >= m_nYGridSize)
      nY0 = m_nYGridSize - 1;

   else if (nY0 < 0)
      nY0 = 0;

   // OK let's go
   int const nDiffX = nX0 - nX1;
   int const nDiffY = nY0 - nY1;

   if (nDiffX == 0)
   {
      // The two points have the same x coordinates, so we just need to constrain the y co-ord
      if (nY1 < nY0)
      {
         nY1 = -1;

         do
         {
            nY1++;
         } while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
      else
      {
         nY1 = m_nYGridSize;

         do
         {
            nY1--;
         } while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
   }
   else if (nDiffY == 0)
   {
      // The two points have the same y coordinates, so we just need to constrain the x co-ord
      if (nX1 < nX0)
      {
         nX1 = -1;

         do
         {
            nX1++;
         } while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
      else
      {
         nX1 = m_nXGridSize;

         do
         {
            nX1--;
         } while (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue());

         return;
      }
   }
   else
   {
      // The two points have different x coordinates and different y coordinates, so we have to work harder. First find which of the coordinates is the greatest distance outside the grid, and constrain that co-ord for efficiency (since this will reduce the number of times round the loop). Note that both may be inside the grid, if the incorrect co-ord is in the invalid margin, in which case arbitrarily contrain the x co-ord
      int nXDistanceOutside = 0;
      int nYDistanceOutside = 0;

      if (nX1 < 0)
         nXDistanceOutside = -nX1;
      else if (nX1 >= m_nXGridSize)
         nXDistanceOutside = nX1 - m_nXGridSize + 1;

      if (nY1 < 0)
         nYDistanceOutside = -nY1;
      else if (nY1 >= m_nYGridSize)
         nXDistanceOutside = nY1 - m_nYGridSize + 1;

      if (nXDistanceOutside >= nYDistanceOutside)
      {
         // Constrain the x co-ord
         if (nX1 < nX0)
         {
            // The incorrect x co-ord is less than the correct x co-ord: constrain it and find the y co-ord
            nX1 = -1;

            do
            {
               nX1++;

               nY1 = nY0 + nRound(((nX1 - nX0) * nDiffY) / static_cast<double>(nDiffX));
            } while ((nY1 < 0) || (nY1 >= m_nYGridSize) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }
         else
         {
            // The incorrect x co-ord is greater than the correct x-co-ord: constrain it and find the y co-ord
            nX1 = m_nXGridSize;

            do
            {
               nX1--;

               nY1 = nY0 + nRound(((nX1 - nX0) * nDiffY) / static_cast<double>(nDiffX));
            } while ((nY1 < 0) || (nY1 >= m_nYGridSize) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }
      }
      else
      {
         // Constrain the y co-ord
         if (nY1 < nY0)
         {
            // The incorrect y co-ord is less than the correct y-co-ord: constrain it and find the x co-ord
            nY1 = -1;

            do
            {
               nY1++;

               nX1 = nX0 + nRound(((nY1 - nY0) * nDiffX) / static_cast<double>(nDiffY));
            } while ((nX1 < 0) || (nX1 >= m_nXGridSize) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }
         else
         {
            // The incorrect y co-ord is greater than the correct y co-ord: constrain it and find the x co-ord
            nY1 = m_nYGridSize;

            do
            {
               nY1--;

               nX1 = nX0 +
                     nRound(((nY1 - nY0) * nDiffX) / static_cast<double>(nDiffY));
            } while ((nX1 < 0) || (nX1 >= m_nXGridSize) || (m_pRasterGrid->m_Cell[nX1][nY1].bBasementElevIsMissingValue()));

            return;
         }
      }
   }
}

//===============================================================================================================================
//! Constrains the supplied angle to be within 0 and 360 degrees
//===============================================================================================================================
double CSimulation::dKeepWithin360(double const dAngle)
{
   double dNewAngle = dAngle;

   // Sort out -ve angles
   while (dNewAngle < 0)
      dNewAngle += 360;

   // Sort out angles > 360
   while (dNewAngle > 360)
      dNewAngle -= 360;

   return dNewAngle;
}

//===============================================================================================================================
//! Returns a point (external CRS) which is the average of (i.e. is midway between) two other external CRS points
//===============================================================================================================================
CGeom2DPoint CSimulation::PtAverage(CGeom2DPoint const* pPt1, CGeom2DPoint const* pPt2)
{
   double const dPt1X = pPt1->dGetX();
   double const dPt1Y = pPt1->dGetY();
   double const dPt2X = pPt2->dGetX();
   double const dPt2Y = pPt2->dGetY();
   double const dPtAvgX = (dPt1X + dPt2X) / 2;
   double const dPtAvgY = (dPt1Y + dPt2Y) / 2;

   return CGeom2DPoint(dPtAvgX, dPtAvgY);
}

// //===============================================================================================================================
// //! Returns an integer point (grid CRS) which is the approximate average of (i.e. is midway between) two other grid CRS integer points
// //===============================================================================================================================
// CGeom2DIPoint CSimulation::PtiAverage(CGeom2DIPoint const* pPti1, CGeom2DIPoint const* pPti2)
// {
//    int const nPti1X = pPti1->nGetX();
//    int const nPti1Y = pPti1->nGetY();
//    int const nPti2X = pPti2->nGetX();
//    int const nPti2Y = pPti2->nGetY();
//    int const nPtiAvgX = (nPti1X + nPti2X) / 2;
//    int const nPtiAvgY = (nPti1Y + nPti2Y) / 2;
//
//    return CGeom2DIPoint(nPtiAvgX, nPtiAvgY);
// }

//===============================================================================================================================
//! Returns an integer point (grid CRS) which is the weighted average of two other grid CRS integer points. The weight must be <= 1, if the weight is < 0.5 then the output point is closer to the first point, if the weight is > 0.5 then the output point is closer to the second point
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiWeightedAverage(CGeom2DIPoint const* pPti1, CGeom2DIPoint const* pPti2, double const dWeight)
{
   int const nPti1X = pPti1->nGetX();
   int const nPti1Y = pPti1->nGetY();
   int const nPti2X = pPti2->nGetX();
   int const nPti2Y = pPti2->nGetY();
   double const dOtherWeight = 1.0 - dWeight;

   int const nPtiWeightAvgX = nRound((dWeight * nPti2X) + (dOtherWeight * nPti1X));
   int const nPtiWeightAvgY = nRound((dWeight * nPti2Y) + (dOtherWeight * nPti1Y));

   return CGeom2DIPoint(nPtiWeightAvgX, nPtiWeightAvgY);
}

//===============================================================================================================================
//! Returns a point (external CRS) which is the average of a vector of external CRS points
//===============================================================================================================================
CGeom2DPoint CSimulation::PtAverage(vector<CGeom2DPoint>* pVIn)
{
   int const nSize = static_cast<int>(pVIn->size());

   if (nSize == 0)
      return CGeom2DPoint(DBL_NODATA, DBL_NODATA);

   double dAvgX = 0;
   double dAvgY = 0;

   for (int n = 0; n < nSize; n++)
   {
      dAvgX += pVIn->at(n).dGetX();
      dAvgY += pVIn->at(n).dGetY();
   }

   dAvgX /= nSize;
   dAvgY /= nSize;

   return CGeom2DPoint(dAvgX, dAvgY);
}

// //===============================================================================================================================
// //! Returns a point (grid CRS) which is the average of a vector of grid CRS points
// //===============================================================================================================================
// CGeom2DIPoint CSimulation::PtiAverage(vector<CGeom2DIPoint>* pVIn)
// {
// int nSize = static_cast<int>(pVIn->size());
// if (nSize == 0)
// return CGeom2DIPoint(INT_NODATA, INT_NODATA);
//
// double dAvgX = 0;
// double dAvgY = 0;
//
// for (int n = 0; n < nSize; n++)
// {
// dAvgX += pVIn->at(n).nGetX();
// dAvgY += pVIn->at(n).nGetY();
// }
//
// dAvgX /= nSize;
// dAvgY /= nSize;
//
// return CGeom2DIPoint(nRound(dAvgX), nRound(dAvgY));
// }

//===============================================================================================================================
//! Returns an integer point (grid CRS) which is the centroid of a polygon, given by a vector of grid CRS points. From https://stackoverflow.com/questions/2792443/finding-the-centroid-of-a-polygon
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiPolygonCentroid(vector<CGeom2DIPoint>* pVIn)
{
   CGeom2DIPoint PtiCentroid(0, 0);
   int const nSize = static_cast<int>(pVIn->size());
   int nX0 = 0; // Current vertex X
   int nY0 = 0; // Current vertex Y
   int nX1 = 0; // Next vertex X
   int nY1 = 0; // Next vertex Y

   double dA = 0; // Partial signed area
   double dSignedArea = 0.0;

   // For all vertices except last
   for (int i = 0; i < nSize - 1; ++i)
   {
      nX0 = pVIn->at(i).nGetX();
      nY0 = pVIn->at(i).nGetY();
      nX1 = pVIn->at(i + 1).nGetX();
      nY1 = pVIn->at(i + 1).nGetY();

      dA = (nX0 * nY1) - (nX1 * nY0);
      dSignedArea += dA;
      PtiCentroid.AddXAddY((nX0 + nX1) * dA, (nY0 + nY1) * dA);
   }

   // Do last vertex separately to avoid performing an expensive modulus operation in each iteration
   nX0 = pVIn->at(nSize - 1).nGetX();
   nY0 = pVIn->at(nSize - 1).nGetY();
   nX1 = pVIn->at(0).nGetX();
   nY1 = pVIn->at(0).nGetY();

   dA = (nX0 * nY1) - (nX1 * nY0);
   dSignedArea += dA;
   PtiCentroid.AddXAddY((nX0 + nX1) * dA, (nY0 + nY1) * dA);

   dSignedArea *= 0.5;
   PtiCentroid.DivXDivY(6.0 * dSignedArea, 6.0 * dSignedArea);

   return PtiCentroid;
}

//===============================================================================================================================
//   Returns a vector which is perpendicular to an existing vector
//===============================================================================================================================
// vector<CGeom2DPoint> CSimulation::VGetPerpendicular(CGeom2DPoint const*
// PtStart, CGeom2DPoint const* PtNext, double const dDesiredLength, int const
// nHandedness)
// {
//    // Returns a two-point vector which passes through PtStart with a scaled
//    length
// double dXLen = PtNext->dGetX() - PtStart->dGetX();
// double dYLen = PtNext->dGetY() - PtStart->dGetY();
//
// double dLength = hypot(dXLen, dYLen);
// double dScaleFactor = dDesiredLength / dLength;
//
//    // The difference vector is (dXLen, dYLen), so the perpendicular
//    difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
// CGeom2DPoint EndPt;
// if (nHandedness == RIGHT_HANDED)
// {
// EndPt.SetX(PtStart->dGetX() + (dScaleFactor * dYLen));
// EndPt.SetY(PtStart->dGetY() - (dScaleFactor * dXLen));
// }
// else
// {
// EndPt.SetX(PtStart->dGetX() - (dScaleFactor * dYLen));
// EndPt.SetY(PtStart->dGetY() + (dScaleFactor * dXLen));
// }
//
// vector<CGeom2DPoint> VNew;
// VNew.push_back(*PtStart);
// VNew.push_back(EndPt);
// return VNew;
// }

// //===============================================================================================================================
// //! Returns a CGeom2DPoint which is the 'other' point of a two-point vector passing through PtStart, and which is perpendicular to the two-point vector from PtStart to PtNext
// //===============================================================================================================================
// CGeom2DPoint CSimulation::PtGetPerpendicular(CGeom2DPoint const* PtStart, CGeom2DPoint const* PtNext, double const dDesiredLength, int const nHandedness)
// {
//    double const dXLen = PtNext->dGetX() - PtStart->dGetX();
//    double const dYLen = PtNext->dGetY() - PtStart->dGetY();
//    double dLength;
//
//    if (bFPIsEqual(dXLen, 0.0, TOLERANCE))
//       dLength = dYLen;
//    else if (bFPIsEqual(dYLen, 0.0, TOLERANCE))
//       dLength = dXLen;
//    else
//       dLength = hypot(dXLen, dYLen);
//
//    double const dScaleFactor = dDesiredLength / dLength;
//
//    // The difference vector is (dXLen, dYLen), so the perpendicular difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
//    CGeom2DPoint EndPt;
//
//    if (nHandedness == RIGHT_HANDED)
//    {
//       EndPt.SetX(PtStart->dGetX() + (dScaleFactor * dYLen));
//       EndPt.SetY(PtStart->dGetY() - (dScaleFactor * dXLen));
//    }
//    else
//    {
//       EndPt.SetX(PtStart->dGetX() - (dScaleFactor * dYLen));
//       EndPt.SetY(PtStart->dGetY() + (dScaleFactor * dXLen));
//    }
//
//    return EndPt;
// }

//===============================================================================================================================
//! Returns a CGeom2DIPoint (grid CRS) which is the inland (i.e. on the right of a left-handed coast) point of a two-point vector passing through PtiStart, and which is perpendicular to the two-point vector from PtiStart to PtiNext
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiGetPerpendicular(CGeom2DIPoint const* PtiStart, CGeom2DIPoint const* PtiNext, double const dDesiredLength, int const nHandedness)
{
   double const dXLen = PtiNext->nGetX() - PtiStart->nGetX();
   double const dYLen = PtiNext->nGetY() - PtiStart->nGetY();
   double dLength;

   if (bFPIsEqual(dXLen, 0.0, TOLERANCE))
      dLength = dYLen;
   else if (bFPIsEqual(dYLen, 0.0, TOLERANCE))
      dLength = dXLen;
   else
      dLength = hypot(dXLen, dYLen);

   double const dScaleFactor = dDesiredLength / dLength;

   // The difference vector is (dXLen, dYLen), so the perpendicular difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
   CGeom2DIPoint EndPti;

   // Choose the inland point
   if (nHandedness == RIGHT_HANDED)
   {
      EndPti.SetX(PtiStart->nGetX() + nRound(dScaleFactor * dYLen));
      EndPti.SetY(PtiStart->nGetY() - nRound(dScaleFactor * dXLen));
   }
   else
   {
      EndPti.SetX(PtiStart->nGetX() - nRound(dScaleFactor * dYLen));
      EndPti.SetY(PtiStart->nGetY() + nRound(dScaleFactor * dXLen));
   }

   return EndPti;
}

//===============================================================================================================================
//! Returns a CGeom2DIPoint (grid CRS) which is which is the inland (i.e. on the right of a left-handed coast) point of a two-point vector passing through [nStartX][nStartY], and which is perpendicular to the two-point vector from [nStartX][nStartY] to [nNextX][nNextY]
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiGetPerpendicular(int const nStartX, int const nStartY, int const nNextX, int const nNextY, double const dDesiredLength, int const nHandedness)
{
   double const dXLen = nNextX - nStartX;
   double const dYLen = nNextY - nStartY;
   double dLength;

   if (bFPIsEqual(dXLen, 0.0, TOLERANCE))
      dLength = dYLen;
   else if (bFPIsEqual(dYLen, 0.0, TOLERANCE))
      dLength = dXLen;
   else
      dLength = hypot(dXLen, dYLen);

   double const dScaleFactor = dDesiredLength / dLength;

   // The difference vector is (dXLen, dYLen), so the perpendicular difference vector is (-dYLen, dXLen) or (dYLen, -dXLen)
   CGeom2DIPoint EndPti;

   // Chhose the inland point
   if (nHandedness == RIGHT_HANDED)
   {
      EndPti.SetX(nStartX + nRound(dScaleFactor * dYLen));
      EndPti.SetY(nStartY - nRound(dScaleFactor * dXLen));
   }
   else
   {
      EndPti.SetX(nStartX - nRound(dScaleFactor * dYLen));
      EndPti.SetY(nStartY + nRound(dScaleFactor * dXLen));
   }

   return EndPti;
}

//===============================================================================================================================
//! Returns the signed angle BAC (in radians) subtended between three CGeom2DIPoints B A C. From http://stackoverflow.com/questions/3057448/angle-between-3-vertices
//===============================================================================================================================
double CSimulation::dAngleSubtended(CGeom2DIPoint const* pPtiA, CGeom2DIPoint const* pPtiB, CGeom2DIPoint const* pPtiC)
{
   double const dXDistBtoA = pPtiB->nGetX() - pPtiA->nGetX();
   double const dYDistBtoA = pPtiB->nGetY() - pPtiA->nGetY();
   double const dXDistCtoA = pPtiC->nGetX() - pPtiA->nGetX();
   double const dYDistCtoA = pPtiC->nGetY() - pPtiA->nGetY();
   double const dDotProduct = dXDistBtoA * dXDistCtoA + dYDistBtoA * dYDistCtoA;
   double const dPseudoCrossProduct = dXDistBtoA * dYDistCtoA - dYDistBtoA * dXDistCtoA;
   double const dAngle = atan2(dPseudoCrossProduct, dDotProduct);

   return dAngle;
}

//===============================================================================================================================
//! Checks whether the selected raster GDAL driver supports file creation, 32-bit doubles, etc.
//===============================================================================================================================
bool CSimulation::bCheckRasterGISOutputFormat(void)
{
   // Register all available GDAL raster and vector drivers (GDAL 2)
   GDALAllRegister();

   // If the user hasn't specified a GIS output format, assume that we will use the same GIS format as the input basement DEM
   if (m_strRasterGISOutFormat.empty())
      m_strRasterGISOutFormat = m_strGDALBasementDEMDriverCode;

   // Load the raster GDAL driver
   GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName(m_strRasterGISOutFormat.c_str());

   if (NULL == pDriver)
   {
      // Can't load raster GDAL driver. Incorrectly specified?
      cerr << ERR << "Unknown raster GIS output format '" << m_strRasterGISOutFormat << "'." << endl;
      return false;
   }

   // Get the metadata for this raster driver
   char **papszMetadata = pDriver->GetMetadata();

   // for (int i = 0; papszMetadata[i] != NULL; i++)
   // cout << papszMetadata[i] << endl;
   // cout << endl;

   // Need to test if this is a raster driver
   if (!CSLFetchBoolean(papszMetadata, GDAL_DCAP_RASTER, false))
   {
      // This is not a raster driver
      cerr << ERR << "GDAL driver '" << m_strRasterGISOutFormat << "' is not a raster driver. Choose another format." << endl;
      return false;
   }

   // This driver is OK, so store its longname and the default file extension
   string strTmp = CSLFetchNameValue(papszMetadata, "DMD_LONGNAME");
   m_strGDALRasterOutputDriverLongname = strTrim(&strTmp);
   strTmp = CSLFetchNameValue(papszMetadata, "DMD_EXTENSIONS");      // Note DMD_EXTENSION (no S, is a single value) appears not to be implemented for newer drivers
   strTmp = strTrim(&strTmp);

   // We have a space-separated list of one or more file extensions: use the first extension in the list
   long unsigned int const nPos = strTmp.find(SPACE);

   if (nPos == string::npos)
   {
      // No space i.e. just one extension
      m_strGDALRasterOutputDriverExtension = strTmp;
   }

   else
   {
      // There's a space, so we must have more than one extension
      m_strGDALRasterOutputDriverExtension = strTmp.substr(0, nPos);
   }

   // Set up any defaults for raster files that are created using this driver
   SetRasterFileCreationDefaults();

   // Now do various tests of the driver's capabilities
   if (!CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, false))
   {
      // This raster driver does not support the Create() method, does it support CreateCopy()?
      if (!CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATECOPY, false))
      {
         cerr << ERR << "Cannot write using raster GDAL driver '" << m_strRasterGISOutFormat << " since neither Create() or CreateCopy() are supported'. Choose another GDAL raster format." << endl;
         return false;
      }

      // Can't use Create() but can use CreateCopy()
      m_bGDALCanCreate = false;
   }

   // Next, test to see what data types the driver can write and from this, work out the largest int and float we can write
   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Float"))
   {
      m_bGDALCanWriteFloat = true;
      m_GDALWriteFloatDataType = GDT_Float32;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "UInt32"))
   {
      m_bGDALCanWriteInt32 = true;

      m_GDALWriteIntDataType = GDT_UInt32;
      m_lGDALMaxCanWrite = UINT32_MAX;
      m_lGDALMinCanWrite = 0;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_UInt32;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Int32"))
   {
      m_bGDALCanWriteInt32 = true;

      m_GDALWriteIntDataType = GDT_Int32;
      m_lGDALMaxCanWrite = INT32_MAX;
      m_lGDALMinCanWrite = INT32_MIN;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_Int32;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "UInt16"))
   {
      m_bGDALCanWriteInt32 = false;

      m_GDALWriteIntDataType = GDT_UInt16;
      m_lGDALMaxCanWrite = UINT16_MAX;
      m_lGDALMinCanWrite = 0;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_UInt16;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Int16"))
   {
      m_bGDALCanWriteInt32 = false;

      m_GDALWriteIntDataType = GDT_Int16;
      m_lGDALMaxCanWrite = INT16_MAX;
      m_lGDALMinCanWrite = INT16_MIN;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_Int16;

      return true;
   }

   if (strstr(CSLFetchNameValue(papszMetadata, "DMD_CREATIONDATATYPES"), "Byte"))
   {
      m_bGDALCanWriteInt32 = false;

      m_GDALWriteIntDataType = GDT_Byte;
      m_lGDALMaxCanWrite = UINT8_MAX;
      m_lGDALMinCanWrite = 0;

      if (! m_bGDALCanWriteFloat)
         m_GDALWriteFloatDataType = GDT_Byte;

      return true;
   }

   // This driver does not even support byte output
   cerr << ERR << "Cannot write using raster GDAL driver '" << m_strRasterGISOutFormat << ", not even byte output is supported'. Choose another GIS raster format." << endl;
   return false;
}

//===============================================================================================================================
//! Checks whether the selected vector GDAL/OGR driver supports file creation etc.
//===============================================================================================================================
bool CSimulation::bCheckVectorGISOutputFormat(void)
{
   // Load the vector GDAL driver (this assumes that GDALAllRegister() has already been called)
   GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName(m_strVectorGISOutFormat.c_str());

   if (NULL == pDriver)
   {
      // Can't load vector GDAL driver. Incorrectly specified?
      cerr << ERR << "Unknown vector GIS output format '" << m_strVectorGISOutFormat << "'." << endl;
      return false;
   }

   // Get the metadata for this vector driver
   char **papszMetadata = pDriver->GetMetadata();

   // For GDAL2, need to test if this is a vector driver
   if (!CSLFetchBoolean(papszMetadata, GDAL_DCAP_VECTOR, false))
   {
      // This is not a vector driver
      cerr << ERR << "GDAL driver '" << m_strVectorGISOutFormat << "' is not a vector driver. Choose another format." << endl;
      return false;
   }

   if (!CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, false))
   {
      // Driver does not support create() method
      cerr << ERR << "Cannot write vector GIS files using GDAL driver '" << m_strRasterGISOutFormat << "'. Choose another format." << endl;
      return false;
   }

   // Driver is OK, now set some options for individual drivers
   if (m_strVectorGISOutFormat == "ESRI Shapefile")
   {
      // Set this, so that just a single dataset-with-one-layer shapefile is created, rather than a directory (see http://www.gdal.org/ogr/drv_shapefile.html)
      m_strOGRVectorOutputExtension = ".shp";
   }

   else if (m_strVectorGISOutFormat == "geojson")
   {
      m_strOGRVectorOutputExtension = ".geojson";
   }

   else if (m_strVectorGISOutFormat == "gpkg")
   {
      m_strOGRVectorOutputExtension = ".gpkg";
   }

   // TODO 033 Others

   return true;
}

//===============================================================================================================================
//! The bSaveAllRasterGISFiles member function saves the raster GIS files using values from the RasterGrid array
//===============================================================================================================================
bool CSimulation::bSaveAllRasterGISFiles(void)
{
   // Increment file number
   m_nGISSave++;

   // Set for next save
   if (m_bSaveRegular)
   {
      m_dRegularSaveTime += m_dRegularSaveInterval;
   }
   else
   {
      if (m_nThisSave < m_nUSave - 1)
      {
         // Still have user-defined save times remaining
         m_nThisSave++;
      }
      else
      {
         // Finished user-defined times, switch to regular interval using last value as interval
         double dLastInterval;

         if (m_nUSave > 1)
            dLastInterval = m_dUSaveTime[m_nUSave - 1] - m_dUSaveTime[m_nUSave - 2];
         else
            dLastInterval = m_dUSaveTime[m_nUSave - 1];

         m_dRegularSaveTime = m_dSimElapsed + dLastInterval;
         m_dRegularSaveInterval = dLastInterval;
         m_bSaveRegular = true;
      }
   }

   if (m_bSedIncTalusTopSurfSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_SED_TOP_INC_TALUS_ELEV, &RASTER_PLOT_SED_TOP_INC_TALUS_ELEV_TITLE))
         return false;

   if (m_bTopSurfIncSeaSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_TOP_ELEV_INC_SEA, &RASTER_PLOT_TOP_ELEV_INC_SEA_TITLE))
         return false;

   if (m_bTalusSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_TALUS, &RASTER_PLOT_TALUS_TITLE))
         return false;

   if (m_bConsSedSlopeSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_CONS_SED_SLOPE, &RASTER_PLOT_CONS_SED_SLOPE_TITLE))
         return false;

   // if (m_bSlopeSaveForCliffToe)
   //    if (! bWriteRasterGISFile(RASTER_PLOT_SLOPE_FOR_CLIFF_TOE, &RASTER_PLOT_SLOPE_FOR_CLIFF_TOE_TITLE))
   //       return false;
   //
   // if (m_bCliffToeSave)
   //    if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_TOE, &RASTER_PLOT_CLIFF_TOE_TITLE))
   //       return false;

   if (m_bSeaDepthSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_SEA_DEPTH, &RASTER_PLOT_SEA_DEPTH_TITLE))
         return false;

   if (m_bWaveHeightSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_WAVE_HEIGHT, &RASTER_PLOT_WAVE_HEIGHT_TITLE))
         return false;

   if (m_bWaveAngleSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_WAVE_ANGLE, &RASTER_PLOT_WAVE_ANGLE_TITLE))
         return false;

   // Don't write platform erosion files if there is no platform erosion
   if (m_bDoShorePlatformErosion)
   {
      if (m_bPotentialPlatformErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_POTENTIAL_PLATFORM_EROSION, &RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_TITLE))
            return false;

      if (m_bActualPlatformErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_ACTUAL_PLATFORM_EROSION, &RASTER_PLOT_ACTUAL_PLATFORM_EROSION_TITLE))
            return false;

      if (m_bTotalPotentialPlatformErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION, &RASTER_PLOT_TOTAL_POTENTIAL_PLATFORM_EROSION_TITLE))
            return false;

      if (m_bTotalActualPlatformErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION, &RASTER_PLOT_TOTAL_ACTUAL_PLATFORM_EROSION_TITLE))
            return false;

      if (m_bPotentialPlatformErosionMaskSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK, &RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK_TITLE))
            return false;

      if (m_bBeachProtectionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_BEACH_PROTECTION, &RASTER_PLOT_BEACH_PROTECTION_TITLE))
            return false;

      if (m_bPotentialBeachErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_POTENTIAL_BEACH_EROSION, &RASTER_PLOT_POTENTIAL_BEACH_EROSION_TITLE))
            return false;

      if (m_bActualBeachErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_ACTUAL_BEACH_EROSION, &RASTER_PLOT_ACTUAL_BEACH_EROSION_TITLE))
            return false;

      if (m_bTotalPotentialBeachErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION, &RASTER_PLOT_TOTAL_POTENTIAL_BEACH_EROSION_TITLE))
            return false;

      if (m_bTotalActualBeachErosionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION, &RASTER_PLOT_TOTAL_ACTUAL_BEACH_EROSION_TITLE))
            return false;

      if (m_bBeachDepositionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_BEACH_DEPOSITION, &RASTER_PLOT_BEACH_DEPOSITION_TITLE))
            return false;

      if (m_bTotalBeachDepositionSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_BEACH_DEPOSITION, &RASTER_PLOT_TOTAL_BEACH_DEPOSITION_TITLE))
            return false;
   }

   if (m_bLandformSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_LANDFORM, &RASTER_PLOT_LANDFORM_TITLE))
         return false;

   if (m_bAvgWaveHeightSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_AVG_WAVE_HEIGHT, &RASTER_PLOT_AVG_WAVE_HEIGHT_TITLE))
         return false;

   if (m_bAvgWaveAngleSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_AVG_WAVE_ANGLE, &RASTER_PLOT_AVG_WAVE_ANGLE_TITLE))
         return false;

   if (m_bAvgSeaDepthSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_AVG_SEA_DEPTH, &RASTER_PLOT_AVG_SEA_DEPTH_TITLE))
         return false;

   if (m_bSedimentInput && m_bSedimentInputEventSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_SEDIMENT_INPUT, &RASTER_PLOT_SEDIMENT_INPUT_EVENT_TITLE))
         return false;

   // Don't write suspended sediment files if there is no fine sediment
   if (m_bHaveFineSediment)
   {
      if (m_bSuspSedSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_SUSPENDED_SEDIMENT, &RASTER_PLOT_SUSPENDED_SEDIMENT_TITLE))
            return false;

      if (m_bAvgSuspSedSave)
         if (! bWriteRasterGISFile(RASTER_PLOT_AVG_SUSPENDED_SEDIMENT, &RASTER_PLOT_AVG_SUSPENDED_SEDIMENT_TITLE))
            return false;
   }

   if (m_bBasementElevSave)
      if (! bWriteRasterGISFile(RASTER_PLOT_BASEMENT_ELEVATION, &RASTER_PLOT_BASEMENT_ELEVATION_TITLE))
         return false;

   for (int nLayer = 0; nLayer < m_nLayers; nLayer++)
   {
      if (m_bHaveFineSediment && m_bFineUnconsSedSave)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT, &RASTER_PLOT_FINE_UNCONSOLIDATED_SEDIMENT_TITLE, nLayer))
            return false;
      }

      if (m_bHaveSandSediment && m_bSandUnconsSedSave)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT, &RASTER_PLOT_SAND_UNCONSOLIDATED_SEDIMENT_TITLE, nLayer))
            return false;
      }

      if (m_bHaveCoarseSediment && m_bCoarseUnconsSedSave)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT, &RASTER_PLOT_COARSE_UNCONSOLIDATED_SEDIMENT_TITLE, nLayer))
            return false;
      }

      if (m_bHaveFineSediment && m_bFineConsSedSave)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT, &RASTER_PLOT_FINE_CONSOLIDATED_SEDIMENT_TITLE, nLayer))
            return false;
      }

      if (m_bHaveSandSediment && m_bSandConsSedSave)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT, &RASTER_PLOT_SAND_CONSOLIDATED_SEDIMENT_TITLE, nLayer))
            return false;
      }

      if (m_bHaveCoarseSediment && m_bCoarseConsSedSave)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT, &RASTER_PLOT_COARSE_CONSOLIDATED_SEDIMENT_TITLE, nLayer))
            return false;
      }
   }

   if (m_bSliceSave)
   {
      for (int i = 0; i < static_cast<int>(m_VdSliceElev.size()); i++)
      {
         if (! bWriteRasterGISFile(RASTER_PLOT_SLICE, &RASTER_PLOT_SLICE_TITLE, 0, m_VdSliceElev[i]))
            return false;
      }
   }

   if (m_bRasterCoastlineSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_COAST, &RASTER_PLOT_COAST_TITLE))
         return false;
   }

   if (m_bRasterNormalProfileSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_NORMAL_PROFILE, &RASTER_PLOT_NORMAL_PROFILE_TITLE))
         return false;
   }

   if (m_bActiveZoneSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_ACTIVE_ZONE, &RASTER_PLOT_ACTIVE_ZONE_TITLE))
         return false;
   }

   // Don't write cliff collapse files if we aren't considering cliff collapse
   if (m_bDoCliffCollapse)
   {
      if (m_bCliffCollapseSave)
      {
         if (m_bHaveFineSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE, &RASTER_PLOT_CLIFF_COLLAPSE_EROSION_FINE_TITLE))
               return false;
         }

         if (m_bHaveSandSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND, &RASTER_PLOT_CLIFF_COLLAPSE_EROSION_SAND_TITLE))
               return false;
         }

         if (m_bHaveCoarseSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE, &RASTER_PLOT_CLIFF_COLLAPSE_EROSION_COARSE_TITLE))
               return false;
         }

         if (m_bCliffNotchAllSave)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_NOTCH_ALL, &RASTER_PLOT_CLIFF_NOTCH_ALL_TITLE))
               return false;
         }

         if (m_bCliffCollapseTimestepSave)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP, &RASTER_PLOT_CLIFF_COLLAPSE_TIMESTEP_TITLE))
               return false;
         }
      }

      if (m_bTotCliffCollapseSave)
      {
         if (m_bHaveFineSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE, &RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_TITLE))
               return false;
         }

         if (m_bHaveSandSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND, &RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_TITLE))
               return false;
         }

         if (m_bHaveCoarseSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE, &RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_TITLE))
               return false;
         }
      }

      if (m_bCliffCollapseDepositionSave)
      {
         if (m_bHaveSandSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND, &RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_SAND_TITLE))
               return false;
         }

         if (m_bHaveCoarseSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE, &RASTER_PLOT_CLIFF_COLLAPSE_DEPOSITION_COARSE_TITLE))
               return false;
         }
      }

      if (m_bTotCliffCollapseDepositionSave)
      {
         if (m_bHaveSandSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND, &RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_TITLE))
               return false;
         }

         if (m_bHaveCoarseSediment)
         {
            if (! bWriteRasterGISFile(RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE, &RASTER_PLOT_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_TITLE))
               return false;
         }
      }
   }

   if (m_bRasterPolygonSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON, &RASTER_PLOT_POLYGON_TITLE))
         return false;
   }

   if (m_bSeaMaskSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_INUNDATION_MASK, &RASTER_PLOT_INUNDATION_MASK_TITLE))
         return false;
   }

   if (m_bBeachMaskSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_BEACH_MASK, &RASTER_PLOT_BEACH_MASK_TITLE))
         return false;
   }

   if (m_bInterventionClassSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_INTERVENTION_CLASS, &RASTER_PLOT_INTERVENTION_CLASS_TITLE))
         return false;
   }

   if (m_bInterventionHeightSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_INTERVENTION_HEIGHT, &RASTER_PLOT_INTERVENTION_HEIGHT_TITLE))
         return false;
   }

   if (m_bShadowZoneCodesSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_SHADOW_ZONE, &RASTER_PLOT_SHADOW_ZONE_TITLE))
         return false;

      if (! bWriteRasterGISFile(RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE, &RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE_TITLE))
         return false;
   }

   if (m_bDeepWaterWaveAngleSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_DEEP_WATER_WAVE_ANGLE, &RASTER_PLOT_DEEP_WATER_WAVE_ANGLE_TITLE))
         return false;
   }

   if (m_bDeepWaterWaveHeightSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT, &RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT_TITLE))
         return false;
   }

   if (m_bDeepWaterWavePeriodSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_DEEP_WATER_WAVE_PERIOD, &RASTER_PLOT_DEEP_WATER_WAVE_PERIOD_TITLE))
         return false;
   }

   if (m_bPolygonUnconsSedUpOrDownDriftSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT, &RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT_TITLE))
         return false;
   }

   if (m_bPolygonUnconsSedGainOrLossSave)
   {
      if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON_GAIN_OR_LOSS, &RASTER_PLOT_POLYGON_GAIN_OR_LOSS_TITLE))
         return false;
   }

   // if (m_bSetupSurgeFloodMaskSave)
   // {
   //    if (! bWriteRasterGISFile(RASTER_PLOT_SETUP_SURGE_FLOOD_MASK, &RASTER_PLOT_SETUP_SURGE_FLOOD_MASK_TITLE))
   //       return false;
   // }

   // if (m_bSetupSurgeRunupFloodMaskSave)
   // {
   //    if (! bWriteRasterGISFile(RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK, &RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK_TITLE))
   //       return false;
   // }
   //
   return true;
}

//===============================================================================================================================
//! The bSaveAllvectorGISFiles member function saves the vector GIS files TODO 081 Choose more files to omit from "usual" vector output
//===============================================================================================================================
bool CSimulation::bSaveAllVectorGISFiles(void)
{
   // Always written
   // if (m_bCoastSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_COAST, &VECTOR_PLOT_COAST_TITLE))
         return false;

      if (! bWriteVectorGISFile(VECTOR_PLOT_COAST_SWL_HIGHEST, &VECTOR_PLOT_COAST_SWL_HIGHEST_TITLE))
         return false;

      if (! bWriteVectorGISFile(VECTOR_PLOT_COAST_SWL_LOWEST, &VECTOR_PLOT_COAST_SWL_LOWEST_TITLE))
         return false;
   }

   if (m_bCliffEdgeSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_CLIFF_EDGE, &VECTOR_PLOT_CLIFF_EDGE_TITLE))
         return false;
   }

   // if (m_bNormalsSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_NORMALS, &VECTOR_PLOT_NORMALS_TITLE))
         return false;
   }

   // if (m_bInvalidNormalsSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_INVALID_NORMALS, &VECTOR_PLOT_INVALID_NORMALS_TITLE))
         return false;
   }

   if (m_bCoastCurvatureSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_COAST_CURVATURE, &VECTOR_PLOT_COAST_CURVATURE_TITLE))
         return false;
   }

   // if (m_bWaveAngleAndHeightSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT, &VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT_TITLE))
         return false;
   }

   if (m_bWaveEnergySinceCollapseSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE, &VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE_TITLE))
         return false;
   }

   if (m_bMeanWaveEnergySave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_MEAN_WAVE_ENERGY, &VECTOR_PLOT_MEAN_WAVE_ENERGY_TITLE))
         return false;
   }

   if (m_bBreakingWaveHeightSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_BREAKING_WAVE_HEIGHT, &VECTOR_PLOT_BREAKING_WAVE_HEIGHT_TITLE))
         return false;
   }

   if (m_bPolygonNodeSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_POLYGON_NODES, &VECTOR_PLOT_POLYGON_NODES_TITLE))
         return false;
   }

   // if (m_bPolygonBoundarySave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_POLYGON_BOUNDARY, &VECTOR_PLOT_POLYGON_BOUNDARY_TITLE))
         return false;
   }

   if (m_bCliffNotchSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_CLIFF_NOTCH_ACTIVE, &VECTOR_PLOT_CLIFF_NOTCH_ACTIVE_TITLE))
         return false;
   }

   if (m_bWaveTransectPointsSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_WAVE_TRANSECT_POINTS, &VECTOR_PLOT_WAVE_TRANSECT_POINTS_TITLE))
         return false;
   }

   // if (m_bShadowBoundarySave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_SHADOW_ZONE_BOUNDARY, &VECTOR_PLOT_SHADOW_ZONE_BOUNDARY_TITLE))
         return false;
   }

   // if (m_bShadowDowndriftBoundarySave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_DOWNDRIFT_ZONE_BOUNDARY, &VECTOR_PLOT_DOWNDRIFT_ZONE_BOUNDARY_TITLE))
         return false;
   }

   if (m_bDeepWaterWaveAngleAndHeightSave)
   {
      if (! bWriteVectorGISFile(VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT, &VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_TITLE))
         return false;
   }

   // if (m_bWaveSetupSave)
   // {
   //    if (! bWriteVectorGISFile(VECTOR_PLOT_WAVE_SETUP, &VECTOR_PLOT_WAVE_SETUP_TITLE))
   //       return false;
   // }
   //
   // if (m_bStormSurgeSave)
   // {
   //    if (! bWriteVectorGISFile(VECTOR_PLOT_STORM_SURGE, &VECTOR_PLOT_STORM_SURGE_TITLE))
   //       return false;
   // }
   //
   // if (m_bRunUpSave)
   // {
   //    if (! bWriteVectorGISFile(VECTOR_PLOT_RUN_UP, &VECTOR_PLOT_RUN_UP_TITLE))
   //       return false;
   // }
   //
   // if (m_bRiverineFlooding && m_bVectorWaveFloodLineSave)
   // {
   //    if (! bWriteVectorGISFile(VECTOR_PLOT_FLOOD_LINE, &VECTOR_PLOT_FLOOD_SWL_SETUP_LINE_TITLE))
   //       return false;
   //
   //    // if (! bWriteVectorGISFile(VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_LINE,
   //    // &VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_LINE_TITLE)) return false;
   //
   //    // if (! bWriteVectorGISFile(VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE,
   //    // &VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_TITLE)) return false;
   // }

   return true;
}

//===============================================================================================================================
//! Finds the max and min values in order to scale raster output if we cannot write doubles
//===============================================================================================================================
void CSimulation::GetRasterOutputMinMax(int const nDataItem, double& dMin, double& dMax, int const nLayer, double const dElev)
{
   // If this is a binary mask layer, we already know the max and min values
   if ((nDataItem == RASTER_PLOT_POTENTIAL_PLATFORM_EROSION_MASK) ||
       (nDataItem == RASTER_PLOT_INUNDATION_MASK) ||
       (nDataItem == RASTER_PLOT_BEACH_MASK) ||
       (nDataItem == RASTER_PLOT_COAST) ||
       (nDataItem == RASTER_PLOT_NORMAL_PROFILE) ||
       (nDataItem == RASTER_PLOT_ACTIVE_ZONE) ||
       (nDataItem == RASTER_PLOT_POLYGON_UPDRIFT_OR_DOWNDRIFT)) /*||
       (nDataItem == RASTER_PLOT_SETUP_SURGE_FLOOD_MASK) ||
       (nDataItem == RASTER_PLOT_SETUP_SURGE_RUNUP_FLOOD_MASK) ||
       (nDataItem == RASTER_PLOT_WAVE_FLOOD_LINE))*/
   {
      dMin = 0;
      dMax = 1;

      return;
   }

   // Not a binary mask layer, so we must find the max and min values
   dMin = DBL_MAX;
   dMax = DBL_MIN;

   double dTmp = 0;

   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         switch (nDataItem)
         {
         case (RASTER_PLOT_SLICE):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dElev);
            break;

         case (RASTER_PLOT_LANDFORM):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory();
            break;

         case (RASTER_PLOT_INTERVENTION_CLASS):
            dTmp = INT_NODATA;

            if (bIsInterventionCell(nX, nY))
               dTmp = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory();

            break;

         case (RASTER_PLOT_INTERVENTION_HEIGHT):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight();
            break;

         case (RASTER_PLOT_POLYGON):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            break;

         case (RASTER_PLOT_BASEMENT_ELEVATION):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBasementElev();
            break;

         case (RASTER_PLOT_SED_TOP_INC_TALUS_ELEV):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus();
            break;

         case (RASTER_PLOT_TALUS):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTalusDepth();
            break;

         case (RASTER_PLOT_TOP_ELEV_INC_SEA):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTopElevIncSea();
            break;

         case (RASTER_PLOT_CONS_SED_SLOPE):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedSlope();
            break;

         case (RASTER_PLOT_SEA_DEPTH):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();
            break;

         case (RASTER_PLOT_AVG_SEA_DEPTH):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSeaDepth() / static_cast<double>(m_ulIter);
            break;

         case (RASTER_PLOT_WAVE_HEIGHT):
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
               dTmp = m_dMissingValue;
            else
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

            break;

         case (RASTER_PLOT_AVG_WAVE_HEIGHT):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveHeight() / static_cast<double>(m_ulIter);
            break;

         case (RASTER_PLOT_WAVE_ANGLE):
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
               dTmp = m_dMissingValue;

            else
               dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();

            break;

         case (RASTER_PLOT_AVG_WAVE_ANGLE):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveAngle() / static_cast<double>(m_ulIter);
            break;

         case (RASTER_PLOT_BEACH_PROTECTION):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor();

            if (! bFPIsEqual(dTmp, DBL_NODATA, TOLERANCE))
               dTmp = 1 - dTmp;     // Output the inverse, seems more intuitive

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
            dTmp = m_pRasterGrid->m_Cell[nX][nY].dGetTotSuspendedSediment() /
                   static_cast<double>(m_ulIter);
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

         case (RASTER_PLOT_SHADOW_ZONE):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();
            break;

         case (RASTER_PLOT_SHADOW_DOWNDRIFT_ZONE):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
            break;

         case (RASTER_PLOT_DEEP_WATER_WAVE_ANGLE):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
            break;

         case (RASTER_PLOT_DEEP_WATER_WAVE_HEIGHT):
            dTmp = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
            break;

         case (RASTER_PLOT_POLYGON_GAIN_OR_LOSS):
            int const nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
            int const nPolyCoast = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonCoastID();

            if (nPoly == INT_NODATA)
               dTmp = m_dMissingValue;
            else
               dTmp = m_VCoast[nPolyCoast].pGetPolygon(nPoly)->dGetBeachDepositionAndSuspensionAllUncons();

            break;
         }

         if (! bFPIsEqual(dTmp, DBL_NODATA, TOLERANCE))
         {
            if (dTmp > dMax)
               dMax = dTmp;

            if (dTmp < dMin)
               dMin = dTmp;
         }
      }
   }
}

//===============================================================================================================================
//! Sets per-driver defaults for raster files created using GDAL
//===============================================================================================================================
void CSimulation::SetRasterFileCreationDefaults(void)
{
   string const strDriver = strToLower(&m_strRasterGISOutFormat);
   string const strComment = "Created by " + PROGRAM_NAME + " for " + PLATFORM + " " + strGetBuild() + " running on " + strGetComputerName();

   // TODO 034 Do these for all commonly-used file types
   if (strDriver == "aaigrid")
   {
   }
   else if (strDriver == "bmp")
   {
   }
   else if (strDriver == "gtiff")
   {
      if (m_bWorldFile)
         m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "TFW", "YES");

      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "NUM_THREADS", "ALL_CPUS");
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMPRESS", "LZW");
   }
   else if (strDriver == "hfa")
   {
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "NBITS", "4");
   }
   else if (strDriver == "jpeg")
   {
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMMENT", strComment.c_str());
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "QUALITY", "95");
   }
   else if (strDriver == "png")
   {
      if (m_bWorldFile)
         m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "WORLDFILE", "YES");

      // m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "TITLE", "This is the title"); m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "DESCRIPTION", "This is a description");
      // m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COPYRIGHT", "This is some copyright statement");
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMMENT", strComment.c_str());
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "NBITS", "4");
   }
   else if (strDriver == "rst")
   {
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMMENT", strComment.c_str());
   }
   else if (strDriver == "geojson")
   {
      // m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "COMMENT", strComment.c_str());
   }
   else if (strDriver == "gpkg")
   {
      // TODO 065 Does GDAL support overwriting raster gpkg files yet?
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "OVERWRITE", "YES");
      m_papszGDALRasterOptions = CSLSetNameValue(m_papszGDALRasterOptions, "USE_TILE_EXTENT", "YES");
   }
   else if (strDriver == "netcdf")
   {
   }
}

//===============================================================================================================================
//! Returns the opposite direction
//===============================================================================================================================
int CSimulation::nGetOppositeDirection(int const nDirection)
{
   switch (nDirection)
   {
   case NORTH:
      return SOUTH;

   case NORTH_EAST:
      return SOUTH_WEST;

   case EAST:
      return WEST;

   case SOUTH_EAST:
      return NORTH_WEST;

   case SOUTH:
      return NORTH;

   case SOUTH_WEST:
      return NORTH_EAST;

   case WEST:
      return EAST;

   case NORTH_WEST:
      return SOUTH_EAST;
   }

   // Should never get here
   return NO_DIRECTION;
}

// //===============================================================================================================================
// //! Given two integer points, calculates the slope and intercept of the line passing through the points
// //===============================================================================================================================
// void CSimulation::GetSlopeAndInterceptFromPoints(CGeom2DIPoint const* pPti1, CGeom2DIPoint const* pPti2, double& dSlope, double& dIntercept)
// {
// int nX1 = pPti1->nGetX();
// int nY1 = pPti1->nGetY();
// int nX2 = pPti2->nGetX();
// int nY2 = pPti2->nGetY();
//
// double dXDiff = nX1 - nX2;
// double dYDiff = nY1 - nY2;
//
// if (bFPIsEqual(dXDiff, 0.0, TOLERANCE))
// dSlope = 0;
// else
// dSlope = dYDiff / dXDiff;
//
// dIntercept = nY1 - (dSlope * nX1);
// }

//===============================================================================================================================
//! Finds the closest point on any coastline to a given point
//===============================================================================================================================
CGeom2DIPoint CSimulation::PtiFindClosestCoastPoint(int const nX, int const nY, int& nCoastFound)
{
   unsigned int nMinSqDist = UINT_MAX;
   CGeom2DIPoint PtiCoastPoint;

   // Do for every coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int j = 0; j < m_VCoast[nCoast].nGetCoastlineSize(); j++)
      {
         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nXCoast = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(j)->nGetX();
         int const nYCoast = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(j)->nGetY();

         // Calculate the squared distance between this point and the given point
         int const nXDist = nX - nXCoast;
         int const nYDist = nY - nYCoast;

         unsigned int const nSqDist = (nXDist * nXDist) + (nYDist * nYDist);

         // Is this the closest so dar?
         if (nSqDist < nMinSqDist)
         {
            nMinSqDist = nSqDist;
            PtiCoastPoint.SetXY(nXCoast, nYCoast);
            nCoastFound = nCoast;
         }
      }
   }

   return PtiCoastPoint;
}

//===============================================================================================================================
//! Finds the number of the closest point on any coastline to a given point, or INT_NODATA in case of error
//===============================================================================================================================
int CSimulation::nFindClosestCoastPoint(int const nX, int const nY, int& nCoastFound)
{
   unsigned int nMinSqDist = UINT_MAX;
   int nCoastPoint = INT_NODATA;

   // Do for every coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int j = 0; j < m_VCoast[nCoast].nGetCoastlineSize(); j++)
      {
         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nXCoast = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(j)->nGetX();
         int const nYCoast = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(j)->nGetY();

         // Calculate the squared distance between this point and the given point
         int const nXDist = nX - nXCoast;
         int const nYDist = nY - nYCoast;

         unsigned int const nSqDist = (nXDist * nXDist) + (nYDist * nYDist);

         // Is this the closest so dar?
         if (nSqDist < nMinSqDist)
         {
            nMinSqDist = nSqDist;
            nCoastPoint = j;
            nCoastFound = nCoast;
         }
      }
   }

   return nCoastPoint;
}

//===============================================================================================================================
//! Given a length in m, this returns the rounded equivalent number of cells
//===============================================================================================================================
int CSimulation::nConvertMetresToNumCells(double const dLen) const
{
   return nRound(dLen / m_dCellSide);
}

//===============================================================================================================================
//! For the straight line between two points A and B, and given another point C, this finds the closest point on the line A-B to the point C. From https://cboard.cprogramming.com/c-programming/155809-find-closest-point-line.html
//===============================================================================================================================
void CSimulation::FindClosestPointOnStraightLine(double const dAx, double const dAy, double const dBx, double const dBy, double const dPx, double const dPy, double& dXRet, double& dYRet)
{
  double const dAPx = dPx - dAx;
  double const dAPy = dPy - dAy;
  double const dABx = dBx - dAx;
  double const dABy = dBy - dAy;
  double const dMagAB2 = (dABx * dABx) + (dABy * dABy);
  double const dABdotAP = (dABx * dAPx) + (dABy * dAPy);
  double const dT = dABdotAP / dMagAB2;

  if ( dT < 0)
  {
      dXRet = dAx;
      dYRet = dAy;
  }
  else if (dT > 1)
  {
      dXRet = dBx;
      dYRet = dBy;
  }
  else
  {
      dXRet = dAx + (dABx * dT);
      dYRet = dAy + (dABy * dT);
  }
}

//===============================================================================================================================
//! Given two different edge cells, returns true if they are adjacent
//===============================================================================================================================
bool CSimulation::bIsAdjacentEdgeCell(CGeom2DIPoint const* pPt1, CGeom2DIPoint const* pPt2)
{
   int const nX1 = pPt1->nGetX();
   int const nY1 = pPt1->nGetY();
   int const nX2 = pPt2->nGetX();
   int const nY2 = pPt2->nGetY();

   if (nX1 == nX2)
   {
      if ((nY1 == nY2 + 1) || (nY1 == nY2 - 1))
         return true;
   }

   if (nY1 == nY2)
   {
      if ((nX1 == nX2 + 1) || (nX1 == nX2 - 1))
         return true;
   }

   return false;
}
