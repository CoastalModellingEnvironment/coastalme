/*!
   \class CGeomILine
   \brief Geometry class used to represent 2D vector integer line objects
   \details TODO 001 This is a more detailed description of the CCGeomLine class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file i_line.h
   \brief Contains CGeomILine definitions
*/

#ifndef ILINE_H
#define ILINE_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "2di_shape.h"

class CGeomILine : public CA2DIShape
{
 private:
 protected:
   void Display(void) override;

 public:
   CGeomILine(void);
   ~CGeomILine(void) override;

   // CGeom2DIPoint* pPtiGetAt(int const);
   // int nGetXAt(int const);
   // int nGetYAt(int const);
   // void SetXAt(int const, int const);
   // void SetYAt(int const, int const);

   bool bIsPresent(int const, int const);
};
#endif // ILINE_H
