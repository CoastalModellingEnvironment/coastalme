/*!
   \file line.cpp
   \brief CGeomLine routines
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
#include <iostream>
using std::ios;

#include "2d_point.h"
#include "line.h"

//! Constructor
CGeomLine::CGeomLine(void)
{
}

//! Overloaded constructor with two points as parameters
CGeomLine::CGeomLine(CGeom2DPoint const* pPt1, CGeom2DPoint const* pPt2)
{
   m_VPoints.push_back(*pPt1);
   m_VPoints.push_back(*pPt2);
}

//! Overloaded constructor with one parameter, this creates a given number of uninitialised points
CGeomLine::CGeomLine(int const nNum)
{
   CGeom2DPoint const pPt;
   m_VPoints.reserve(nNum);

   for (int n = 0; n <= nNum; n++)
      m_VPoints.push_back(pPt);
}

//! Destructor
CGeomLine::~CGeomLine(void)
{
}

//! Returns the X value at a given place in the line
double CGeomLine::dGetXAt(int const n)
{
   return m_VPoints[n].dGetX();
}

//! Returns the Y value at a given place in the line
double CGeomLine::dGetYAt(int const n)
{
   return m_VPoints[n].dGetY();
}

//! Returns the point at a given place in the line
CGeom2DPoint* CGeomLine::pPtGetAt(int const n)
{
   return &m_VPoints[n];
}

// //! Sets the X value at a given place in the line
// void CGeomLine::SetXAt(int const n, double const x)
// {
// m_VPoints[n].SetX(x);
// }

// //! Sets the Y value at a given place in the line
// void CGeomLine::SetYAt(int const n, double const y)
// {
// m_VPoints[n].SetY(y);
// }

// bool CGeomLine::bIsPresent(CGeom2DPoint* Pt)
// {
// if (find(m_VPoints.begin(), m_VPoints.end(), *Pt) != m_VPoints.end())
// return true;
//
// return false;
// }

//! Instantiates the pure virtual function in the abstract parent class, so that CGeomLine is not an abstract class
void CGeomLine::Display(void)
{
}
