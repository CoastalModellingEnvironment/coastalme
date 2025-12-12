/*!
   \class CACoastLandform
   \brief Abstract class, used as a base class for landform objects on the coastline
   \details TODO 001 This is a more detailed description of the CACoastLandform class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file coast_landform.h
   \brief Contains CACoastLandform definitions
*/

#ifndef COASTLANDFORM_H
#define COASTLANDFORM_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/

#include "2di_point.h"
#include "coast.h"

class CRWCoast; // Forward declaration

class CACoastLandform
{
 private:
 protected:
   //! The coast number on which this coast landform sits
   int m_nCoast;

   //! The point on the coast on which this coast landform sits
   int m_nPointOnCoastline;

   //! Landform category code
   int m_nCategory;

   //! Total accumulated wave energy since beginning of simulation
   double m_dTotAccumWaveEnergy;

   //! Pointer to this landform's coast
   CRWCoast* pCoast;

 public:
   CACoastLandform(void);
   virtual ~CACoastLandform(void);

   int nGetCoast(void) const;
   int nGetPointOnCoast(void) const;
   // void SetLandFormCategory(int const);
   int nGetLandFormCategory(void) const;
   CGeom2DIPoint* pPtiGetCellMarkedAsCliff(void) const;
   // void SetTotAccumWaveEnergy(double const);
   void IncTotAccumWaveEnergy(double const);
   double dGetTotAccumWaveEnergy(void) const;

   // Pure virtual function, makes this an abstract class TODO Derived classes need to instantiate this, even tho' their implementation never uses their instantiation and never gets called. Seems a waste of time. Alternative?
   virtual void Display() = 0;
};
#endif // COASTLANDFORM_H
