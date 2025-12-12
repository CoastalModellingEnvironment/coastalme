/*!
   \file 2d_shape.cpp
   \brief Abstract class, used as a base class for 2D objects (line, area, etc.)
   \details Abstract class, used as a base class for 2D objects (line, area, etc.)
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
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
   return m_VPoints[n];
}

// //! Clears this 2D shape
// void CA2DShape::Clear(void)
// {
//    m_VPoints.clear();
// }

//! Resizes the vector which represents this 2D shape
void CA2DShape::Resize(int const nSize)
{
   m_VPoints.resize(nSize);
}

// Returns the number of elements in this 2D shape
int CA2DShape::nGetSize(void) const
{
   return static_cast<int>(m_VPoints.size());
}

// void CA2DShape::InsertAtFront(double const dX, double const dY)
// {
// m_VPoints.insert(m_VPoints.begin(), CGeom2DPoint(dX, dY));
// }

//! Appends a point to this 2D shape
void CA2DShape::Append(CGeom2DPoint const* pPtNew)
{
   m_VPoints.push_back(*pPtNew);
}

//! Appends a point to this 2D shape
void CA2DShape::Append(double const dX, double const dY)
{
   m_VPoints.push_back(CGeom2DPoint(dX, dY));
}

//! Appends a point to this 2D shape only if the point is not the same as the previous point in the vector
void CA2DShape::AppendIfNotPrevious(double const dX, double const dY)
{
   CGeom2DPoint const PtIn(dX, dY);

   if (m_VPoints.empty())
      m_VPoints.push_back(PtIn);

   else if (m_VPoints.back() != &PtIn)
      m_VPoints.push_back(PtIn);
}

//! Returns the last element of this 2D shape
CGeom2DPoint* CA2DShape::pPtBack(void)
{
   return &m_VPoints.back();
}

// void CA2DShape::SetPoints(const vector<CGeom2DPoint>* VNewPoints)
// {
// m_VPoints = *VNewPoints;
// }

// int CA2DShape::nLookUp(CGeom2DPoint* Pt)
// {
// auto it = find(m_VPoints.begin(), m_VPoints.end(), *Pt);
// if (it != m_VPoints.end())
// return it - m_VPoints.begin();
// else
// return -1;
// }

// double CA2DShape::dGetLength(void) const
// {
// int nSize = m_VPoints.size();
//
// if (nSize < 2)
// return -1;
//
// double dLength = 0;
// for (int n = 1; n < nSize; n++)
// {
// double dXlen = m_VPoints[n].dGetX() - m_VPoints[n-1].dGetX();
// double dYlen = m_VPoints[n].dGetY() - m_VPoints[n-1].dGetY();
//
// dLength += hypot(dXlen, dYlen);
// }
//
// return dLength;
// }

//! Returns the address of the vector which represents this 2D shape
vector<CGeom2DPoint>* CA2DShape::pPtVGetPoints(void)
{
   return &m_VPoints;
}

// //! Computes the centroid of this 2D polygon (which may be outside, if this is a concave polygon). From http://stackoverflow.com/questions/2792443/finding-the-centroid-of-a-polygon
// CGeom2DPoint CA2DShape::PtGetCentroid(void)
// {
// int nVertexCount = static_cast<int>(m_VPoints.size());
// double dSignedArea = 0;
// double dCentroidX = 0;
// double dCentroidY = 0;
//
//    // For all vertices
// for (int i = 0; i < nVertexCount; ++i)
// {
// double dXThis = m_VPoints[i].dGetX();
// double dYThis = m_VPoints[i].dGetY();
// double dXNext = m_VPoints[(i+1) % nVertexCount].dGetX();
// double dYNext = m_VPoints[(i+1) % nVertexCount].dGetY();
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
   reverse(m_VPoints.begin(), m_VPoints.end());
}
