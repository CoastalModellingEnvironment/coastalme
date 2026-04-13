/*!
   \file 2d_shape.cpp
   \brief Abstract class, used as a base class for 2D objects (line, area, etc.)
   \details Abstract class, used as a base class for 2D objects (line, area, etc.)
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "2d_point.h"
#include "2d_shape.h"

//! Constructor
CA2DShape::CA2DShape(void)
{
}

//! Destructor
CA2DShape::~CA2DShape(void)
{
}

//! Operator to return one point of this 2D shape
CGeom2DPoint& CA2DShape::operator[](int const n)
{
   // TODO 055 Maybe add a safety check?
   return m_VPtPoints[n];
}

// //! Clears this 2D shape
// void CA2DShape::Clear(void)
// {
//    m_VPtPoints.clear();
// }

//! Resizes the vector which represents this 2D shape
void CA2DShape::Resize(int const nSize)
{
   m_VPtPoints.resize(nSize);
}

// Returns the number of elements in this 2D shape
int CA2DShape::nGetSize(void) const
{
   return static_cast<int>(m_VPtPoints.size());
}

// void CA2DShape::InsertAtFront(double const dX, double const dY)
// {
// m_VPtPoints.insert(m_VPtPoints.begin(), CGeom2DPoint(dX, dY));
// }

//! Appends a point to this 2D shape
void CA2DShape::Append(CGeom2DPoint const* pPtNew)
{
   m_VPtPoints.push_back(*pPtNew);
}

//! Appends a point to this 2D shape
void CA2DShape::Append(double const dX, double const dY)
{
   m_VPtPoints.push_back(CGeom2DPoint(dX, dY));
}

//! Appends a point to this 2D shape only if the point is not the same as the previous point in the vector
void CA2DShape::AppendIfNotPrevious(double const dX, double const dY)
{
   CGeom2DPoint const PtIn(dX, dY);

   if (m_VPtPoints.empty())
      m_VPtPoints.push_back(PtIn);

   else if (m_VPtPoints.back() != &PtIn)
      m_VPtPoints.push_back(PtIn);
}

//! Returns the last element of this 2D shape
CGeom2DPoint* CA2DShape::pPtBack(void)
{
   return &m_VPtPoints.back();
}

// void CA2DShape::SetPoints(const vector<CGeom2DPoint>* VNewPoints)
// {
// m_VPtPoints = *VNewPoints;
// }

bool CA2DShape::bIsPresent(CGeom2DPoint const* pPt)
{
   double const dPtX = pPt->dGetX();
   double const dPtY = pPt->dGetY();

   for (size_t n = 0; n < m_VPtPoints.size(); n++)
   {
      if (bFPIsEqual(dPtX, m_VPtPoints[n].dGetX(), TOLERANCE) && bFPIsEqual(dPtY, m_VPtPoints[n].dGetY(), TOLERANCE))
         return true;
   }

   return false;
}

// double CA2DShape::dGetLength(void) const
// {
// int nSize = m_VPtPoints.size();
//
// if (nSize < 2)
// return -1;
//
// double dLength = 0;
// for (int n = 1; n < nSize; n++)
// {
// double dXlen = m_VPtPoints[n].dGetX() - m_VPtPoints[n-1].dGetX();
// double dYlen = m_VPtPoints[n].dGetY() - m_VPtPoints[n-1].dGetY();
//
// dLength += hypot(dXlen, dYlen);
// }
//
// return dLength;
// }

//! Returns the address of the vector which represents this 2D shape
vector<CGeom2DPoint>* CA2DShape::pPtVGetPoints(void)
{
   return &m_VPtPoints;
}

// //! Computes the centroid of this 2D polygon (which may be outside, if this is a concave polygon). From http://stackoverflow.com/questions/2792443/finding-the-centroid-of-a-polygon
// CGeom2DPoint CA2DShape::PtGetCentroid(void)
// {
// int nVertexCount = static_cast<int>(m_VPtPoints.size());
// double dSignedArea = 0;
// double dCentroidX = 0;
// double dCentroidY = 0;
//
//    // For all vertices
// for (int i = 0; i < nVertexCount; ++i)
// {
// double dXThis = m_VPtPoints[i].dGetX();
// double dYThis = m_VPtPoints[i].dGetY();
// double dXNext = m_VPtPoints[(i+1) % nVertexCount].dGetX();
// double dYNext = m_VPtPoints[(i+1) % nVertexCount].dGetY();
//
// double dA = (dXThis * dYNext) - (dXNext * dYThis);
// dSignedArea += dA;
//
// dCentroidX += (dXThis + dXNext) * dA;
// dCentroidY += (dYThis + dYNext) * dA;
// }
//
// dSignedArea *= 0.5;
// dCentroidX /= (6 * dSignedArea);
// dCentroidY /= (6 * dSignedArea);
//
// return (CGeom2DPoint(dCentroidX, dCentroidY));
// }

//! Reverses the sequence of points in the vector which represents this 2D polygon
void CA2DShape::Reverse(void)
{
   reverse(m_VPtPoints.begin(), m_VPtPoints.end());
}
