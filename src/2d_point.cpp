/*!
   \file 2d_point.cpp
   \brief Geometry class used to represent 2D point objects with floating-point coordinates
   \details The CGeom2DPoint class is used to represent 2D points where the x and y coordinates are floating-point values
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
#include "cme.h"
#include "2d_point.h"

//! Constructor with no parameters (the X and Y coordinates of the new CGeom2DPoint object are set to zero)
CGeom2DPoint::CGeom2DPoint(void)
    : dX(0),
      dY(0)
{
}

//! Constructor with one CGeom2DPoint parameter (for the X and Y coordinates of the new CGeom2DPoint object)
CGeom2DPoint::CGeom2DPoint(CGeom2DPoint const* pPt)
    : dX(pPt->dGetX()),
      dY(pPt->dGetY())
{
}

//! Constructor with two double parameters, for the X and Y coordinates of the new CGeom2DPoint object
CGeom2DPoint::CGeom2DPoint(double const dNewX, double const dNewY)
    : dX(dNewX),
      dY(dNewY)
{
}

//! Returns the CGeom2DPoint object's double X coordinate
double CGeom2DPoint::dGetX(void) const
{
   return dX;
}

//! Returns the CGeom2DPoint object's double Y coordinate
double CGeom2DPoint::dGetY(void) const
{
   return dY;
}

//! The double parameter sets a value for the CGeom2DIPoint object's X coordinate
void CGeom2DPoint::SetX(double const dNewX)
{
   dX = dNewX;
}

//! The double parameter sets a value for the CGeom2DIPoint object's Y coordinate
void CGeom2DPoint::SetY(double const dNewY)
{
   dY = dNewY;
}

// void CGeom2DPoint::SetXY(double const dNewX, double const dNewY)
// {
// dX = dNewX;
// dY = dNewY;
// }

// void CGeom2DPoint::SetXY(CGeom2DPoint const* Pt)
// {
// dX = Pt->dGetX();
// dY = Pt->dGetY();
// }

//! Sets one CGeom2DPoint object equal to another
CGeom2DPoint& CGeom2DPoint::operator=(CGeom2DPoint const* pPt)
{
   dX = pPt->dGetX();
   dY = pPt->dGetY();
   return *this;
}

//! Compares two CGeom2DPoint pointed-to objects for equality
bool CGeom2DPoint::operator==(CGeom2DPoint const* pPt) const
{
   if ((bFPIsEqual(pPt->dGetX(), dX, TOLERANCE)) && (bFPIsEqual(pPt->dGetY(), dY, TOLERANCE)))
      return true;

   return false;
}

//! Compares two CGeom2DPoint objects for equality
bool CGeom2DPoint::operator==(CGeom2DPoint Pt) const
{
   if ((bFPIsEqual(Pt.dGetX(), dX, TOLERANCE)) && (bFPIsEqual(Pt.dGetY(), dY, TOLERANCE)))
      return true;

   return false;
}

//! Compares two CGeom2DPoint pointed-to objects for inequality
bool CGeom2DPoint::operator!=(CGeom2DPoint const* pPt) const
{
   if ((! bFPIsEqual(pPt->dGetX(), dX, TOLERANCE)) || (! bFPIsEqual(pPt->dGetY(), dY, TOLERANCE)))
      return true;

   return false;
}

//! Compares two CGeom2DPoint objects for inequality
bool CGeom2DPoint::operator!=(CGeom2DPoint Pt) const
{
   if ((! bFPIsEqual(Pt.dGetX(), dX, TOLERANCE)) || (! bFPIsEqual(Pt.dGetY(), dY, TOLERANCE)))
      return true;

   return false;
}
