/*!
   \file 2di_shape.cpp
   \brief Abstract class, used as a base class for integer 2D objects (line, area, etc.)
   \details TODO 001 A more detailed description of these routines.
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
   return m_VPoints[n];
}

//! Returns the last integer point from the vector which represents this 2D shape
CGeom2DIPoint& CA2DIShape::Back(void)
{
   return m_VPoints.back();
}

//! Returns the address of the vector which represents this 2D shape
vector<CGeom2DIPoint>* CA2DIShape::pPtiVGetPoints(void)
{
   return &m_VPoints;
}

//! Clears the vector which represents this 2D shape
void CA2DIShape::Clear(void)
{
   m_VPoints.clear();
}

//! Resizes the vector which represents this 2D shape
void CA2DIShape::Resize(const int nSize)
{
   m_VPoints.resize(nSize);
}

//! Returns the number of integer point in the vector which represents this 2D shape
int CA2DIShape::nGetSize(void) const
{
   return static_cast<int>(m_VPoints.size());
}

// void CA2DIShape::InsertAtFront(int const nX, int const nY)
// {
// m_VPoints.insert(m_VPoints.begin(), CGeom2DIPoint(nX, nY));
// }

//! Appends a new integer point to the vector which represents this 2D shape
void CA2DIShape::Append(CGeom2DIPoint const* pPtiNew)
{
   m_VPoints.push_back(*pPtiNew);
}

//! Appends a new integer point to the vector which represents this 2D shape
void CA2DIShape::Append(int const nX, int const nY)
{
   m_VPoints.push_back(CGeom2DIPoint(nX, nY));
}

//! Appends a new integer point to the vector which represents this 2D shape, but only if the point is not the same as the previous point in the vector
void CA2DIShape::AppendIfNotPrevious(int const nX, int const nY)
{
   CGeom2DIPoint const PtiIn(nX, nY);

   if (m_VPoints.empty())
      m_VPoints.push_back(PtiIn);

   else if (m_VPoints.back() != &PtiIn)
      m_VPoints.push_back(PtiIn);
}

//! Appends a new integer point to the vector which represents this 2D shape, but only if the point is not the same as the previous point in the vector
void CA2DIShape::AppendIfNotPrevious(CGeom2DIPoint const* pPtiIn)
{
   if (m_VPoints.empty())
      m_VPoints.push_back(*pPtiIn);

   else if (m_VPoints.back() != pPtiIn)
      m_VPoints.push_back(*pPtiIn);
}

// void CA2DIShape::SetPoints(const vector<CGeom2DIPoint>* VNewPoints)
// {
// m_VPoints = *VNewPoints;
// }

// int CA2DIShape::nLookUp(CGeom2DIPoint* Pti)
// {
// auto it = find(m_VPoints.begin(), m_VPoints.end(), *Pti);
// if (it != m_VPoints.end())
// return it - m_VPoints.begin();
// else
// return -1;
// }
