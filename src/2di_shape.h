/*!
   \class CA2DIShape
   \brief Abstract class, used as a base class for integer 2D objects (line, area, etc.)
   \details TODO 001 This is a more detailed description of the C2IDShape abstract class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file 2di_shape.h
   \brief Contains CA2DIShape definitions
*/

#ifndef C2DISHAPE_H
#define C2DISHAPE_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <vector>
using std::vector;

#include "2di_point.h"

class CA2DIShape
{
 private:
 protected:
   //! The integer points which comprise the integer-coordinate 2D shape
   vector<CGeom2DIPoint> m_VPoints;

   CA2DIShape(void);
   virtual ~CA2DIShape(void);

   void Clear(void);

   virtual void Display() = 0;

 public:
   CGeom2DIPoint& operator[](int const);

   CGeom2DIPoint& Back(void);
   vector<CGeom2DIPoint>* pPtiVGetPoints(void);

   void Resize(const int);
   int nGetSize(void) const;

   // void InsertAtFront(int const, int const);
   void Append(CGeom2DIPoint const*);
   void Append(int const, int const);
   void AppendIfNotPrevious(int const, int const);
   void AppendIfNotPrevious(CGeom2DIPoint const*);

   // void SetPoints(const vector<CGeom2DIPoint>*);
   // int nLookUp(CGeom2DIPoint*);
};
#endif // C2DISHAPE_H
