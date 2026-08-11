/*!
   \file 2di_shape.cpp
   \brief Abstract class, used as a base class for integer 2D objects (line, area, etc.)
   \details TODO 001 A more detailed description of these routines.
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
#include "2di_point.h"
#include "2di_shape.h"

//! Constructor, no parameters
CA2DIShape::CA2DIShape(void)
{
}

//! Destructor
CA2DIShape::~CA2DIShape(void)
{
}

//! Returns one integer point from the vector which represents this 2D shape
CGeom2DIPoint& CA2DIShape::operator[](int const n)
{
   // TODO 055 Maybe add a safety check?
   return m_VPtiPoints[n];
}

//! Returns one integer point from the vector which represents this 2D shape
CGeom2DIPoint& CA2DIShape::at(int const n)
{
   // TODO 055 Maybe add a safety check?
   return m_VPtiPoints[n];

}

//! Returns the last integer point from the vector which represents this 2D shape
CGeom2DIPoint& CA2DIShape::Back(void)
{
   return m_VPtiPoints.back();
}

//! Returns the address of the vector which represents this 2D shape
vector<CGeom2DIPoint>* CA2DIShape::pPtiVGetPoints(void)
{
   return &m_VPtiPoints;
}

//! Clears the vector which represents this 2D shape
void CA2DIShape::Clear(void)
{
   m_VPtiPoints.clear();
}

//! Resizes the vector which represents this 2D shape
void CA2DIShape::Resize(const int nSize)
{
   m_VPtiPoints.resize(nSize);
}

//! Returns the number of integer point in the vector which represents this 2D shape
int CA2DIShape::nGetSize(void) const
{
   return static_cast<int>(m_VPtiPoints.size());
}

// void CA2DIShape::InsertAtFront(int const nX, int const nY)
// {
// m_VPtPoints.insert(m_VPtiPoints.begin(), CGeom2DIPoint(nX, nY));
// }

//! Appends a new integer point to the vector which represents this 2D shape
void CA2DIShape::Append(CGeom2DIPoint const* pPtiNew)
{
   m_VPtiPoints.push_back(*pPtiNew);
}

//! Appends a new integer point to the vector which represents this 2D shape
void CA2DIShape::Append(int const nX, int const nY)
{
   m_VPtiPoints.push_back(CGeom2DIPoint(nX, nY));
}

//! Appends a new integer point to the vector which represents this 2D shape, but only if the point is not the same as the previous point in the vector
void CA2DIShape::AppendIfNotPrevious(int const nX, int const nY)
{
   CGeom2DIPoint const PtiIn(nX, nY);

   if (m_VPtiPoints.empty())
      m_VPtiPoints.push_back(PtiIn);

   else if (m_VPtiPoints.back() != &PtiIn)
      m_VPtiPoints.push_back(PtiIn);
}

//! Appends a new integer point to the vector which represents this 2D shape, but only if the point is not the same as the previous point in the vector
void CA2DIShape::AppendIfNotPrevious(CGeom2DIPoint const* pPtiIn)
{
   if (m_VPtiPoints.empty())
      m_VPtiPoints.push_back(*pPtiIn);

   else if (m_VPtiPoints.back() != pPtiIn)
      m_VPtiPoints.push_back(*pPtiIn);
}

// void CA2DIShape::SetPoints(const vector<CGeom2DIPoint>* VNewPoints)
// {
// m_VPtiPoints = *VNewPoints;
// }

// int CA2DIShape::nLookUp(CGeom2DIPoint* Pti)
// {
// auto it = find(m_VPtiPoints.begin(), m_VPtiPoints.end(), *Pti);
// if (it != m_VPtiPoints.end())
// return it - m_VPtiPoints.begin();
// else
// return -1;
// }
