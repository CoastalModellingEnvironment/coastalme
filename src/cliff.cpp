/*!
   \file cliff.cpp
   \brief CRWCliff routines
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
#include "cliff.h"
#include "coast.h"

//! Constructor with seven parameters and an initialization list
CRWCliff::CRWCliff(CRWCoast* pCoastIn, int const nCoast, int const nPointOnCoast, double const dCellSide, double const dNotchIncisionIn, double const dNotchApexElevIn, double const dAccumWaveEnergyIn)
{
   m_bCliffHasCollapsed = false;

   pCoast = pCoastIn;

   m_nCoast = nCoast;
   m_nPointOnCoastline = nPointOnCoast;
   m_nCategory = LF_CLIFF_ON_COASTLINE;

   m_dMaxNotchIncision = dCellSide;
   m_dNotchIncision = dNotchIncisionIn;
   m_dNotchApexElev = dNotchApexElevIn;
   m_dTotAccumWaveEnergy = dAccumWaveEnergyIn;
}

//! Destructor
CRWCliff::~CRWCliff(void)
{
}

//! Returns the value of the cliff collapse switch
bool CRWCliff::bHasCollapsed(void) const
{
   return m_bCliffHasCollapsed;
}

//! Flags the cliff as having collapsed
void CRWCliff::SetCliffCollapsed(void)
{
   m_bCliffHasCollapsed = true;
}

//! Returns the elevation of the apex of the erosional notch (in external CRS units)
double CRWCliff::dGetNotchApexElev(void) const
{
   return m_dNotchApexElev;
}

//! Sets the elevation of the apex of the erosional notch (in external CRS units)
void CRWCliff::SetNotchApexElev(double const dNewElev)
{
   m_dNotchApexElev = dNewElev;
}

//! Sets the horizontal incision (in external CRS units) of the erosional notch, measured inland from the side of the cell that touches the sea
void CRWCliff::SetNotchIncision(double const dNewIncision)
{
   m_dNotchIncision = dNewIncision;
}

//! Returns the horizontal incision (in external CRS units) of the cliff's erosional notch (the 'overhang')
double CRWCliff::dGetNotchIncision(void) const
{
   return m_dNotchIncision;
}

//! Increases the horizontal incision (in external CRS units) of the erosional notch, measured inland from the side of the cell that touches the sea
void CRWCliff::IncreaseNotchIncision(double const dLenIn)
{
   m_dNotchIncision += dLenIn;

   // Constrain the notch incision, it cannot be greater than the max notch incision
   m_dNotchIncision = tMin(m_dNotchIncision, m_dMaxNotchIncision);

   // assert((m_dMaxNotchIncision - m_dNotchIncision) >=0);
}

//! Returns true if the horizontal incision of the erosional notch exceeds the critical notch incision
bool CRWCliff::bReadyToCollapse(double const dThresholdNotchIncision) const
{
   if (m_dNotchIncision > dThresholdNotchIncision)
      return true;
   else
      return false;
}

//! Instantiates the pure virtual function in the abstract parent class, so that CRWCliff is not an abstract class
void CRWCliff::Display(void)
{
}
