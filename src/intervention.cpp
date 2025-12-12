/*!
   \file intervention.cpp
   \brief CRWIntervention routines
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
#include <assert.h>

#include <iostream>
using std::ios;

#include "cme.h"
#include "coast.h"
#include "intervention.h"

//! Constructor with four parameters
CRWIntervention::CRWIntervention(CRWCoast* pCoastIn, int const nCoast, int const nPointOnCoast, int const nLandCategory)
{
   pCoast = pCoastIn;

   m_nCoast = nCoast;
   m_nPointOnCoastline = nPointOnCoast;
   m_nCategory = nLandCategory;
}

//! Destructor
CRWIntervention::~CRWIntervention(void)
{
}

//! Instantiates the pure virtual function in the abstract parent class, so that CRWIntervention is not an abstract class
void CRWIntervention::Display(void)
{
}
