/*!
   \class CA2DShape
   \brief Abstract class, used as a base class for 2D objects (line, area, etc.)
   \details Abstract class, used as a base class for 2D objects (line, area, etc.)
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file 2d_shape.h
   \brief Contains CA2DShape definitions
*/

#ifndef C2DSHAPE_H
#define C2DSHAPE_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <vector>
using std::vector;

#include <algorithm>
using std::reverse;

#include "2d_point.h"

class CA2DShape
{
 private:
 protected:
   //! The points which comprise the float-coordinate 2D shape
   vector<CGeom2DPoint> m_VPoints;

   CA2DShape(void);
   virtual ~CA2DShape(void);

   // void Clear(void);

   // void InsertAtFront(double const, double const);
   // void SetPoints(const vector<CGeom2DPoint>*);
   // int nLookUp(CGeom2DPoint*);
   // double dGetLength(void) const;
   // CGeom2DPoint PtGetCentroid(void);

   virtual void Display() = 0;

 public:
   void Reverse(void);

   int nGetSize(void) const;
   void Resize(int const);

   void Append(CGeom2DPoint const*);
   void Append(double const, double const);
   void AppendIfNotPrevious(double const, double const);
   CGeom2DPoint* pPtBack(void);

   CGeom2DPoint& operator[](int const);
   vector<CGeom2DPoint>* pPtVGetPoints(void);
};
#endif // C2DSHAPE_H
