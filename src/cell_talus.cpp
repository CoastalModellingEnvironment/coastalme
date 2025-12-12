/*!
   \file cell_talus.cpp
   \brief CRWCellTalus routines
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

#include "cell_talus.h"

//! CRWCellTalus constructor with initialisation list
CRWCellTalus::CRWCellTalus(void)
    : m_dSand(0),
      m_dCoarse(0),
      m_dSandLostThisIter(0),
      m_dCoarseLostThisIter(0),
      m_dSandInputThisIter(0),
      m_dCoarseInputThisIter(0),
      m_dTotSandInput(0),
      m_dTotCoarseInput(0),
      m_dTotSandLost(0),
      m_dTotCoarseLost(0)
{
}

//! Sets this talus object's sand sediment depth equivalent. Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
void CRWCellTalus::SetSandDepth(double const dNewSedDepth)
{
   m_dSand = dNewSedDepth;
}

//! Returns the sand sediment depth equivalent for this talus object
double CRWCellTalus::dGetSandDepth(void) const
{
   return m_dSand;
}

//! Adds sand sediment (depth equivalent) to this talus object object's sand sediment
void CRWCellTalus::AddSandDepth(double const dSedDepthToAdd)
{
   m_dSand += dSedDepthToAdd;
}

//! Sets this talus object object's coarse sediment depth equivalent. Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
void CRWCellTalus::SetCoarseDepth(double const dNewSedDepth)
{
   m_dCoarse = dNewSedDepth;
}

//! Returns the coarse sediment depth equivalent for this talus object object
double CRWCellTalus::dGetCoarseDepth(void) const
{
   return m_dCoarse;
}

//! Adds coarse sediment (depth equivalent) to this talus object object's coarse sediment
void CRWCellTalus::AddCoarseDepth(double const dSedDepthToAdd)
{
   m_dCoarse += dSedDepthToAdd;
}

//! Returns the value for sand talus added during this iteration
double CRWCellTalus::dGetSandAddedThisIter(void) const
{
   return m_dSandInputThisIter;
}

//! Returns the value for sand talus lost during this iteration
double CRWCellTalus::dGetSandLostThisIter(void) const
{
   return m_dSandLostThisIter;
}

//! Returns the value for total sand talus added
double CRWCellTalus::dGetTotSandAdded(void) const
{
   return m_dTotSandInput;
}

//! Returns the value for total sand talus lost
double CRWCellTalus::dGetTotSandLost(void) const
{
   return m_dTotSandLost;
}

//! Returns the value for coarse talus added during this iteration
double CRWCellTalus::dGetCoarseAddedThisIter(void) const
{
   return m_dCoarseInputThisIter;
}

//! Returns the value for coarse talus lost during this iteration
double CRWCellTalus::dGetCoarseLostThisIter(void) const
{
   return m_dCoarseLostThisIter;
}

//! Returns the value for total coarse talus added
double CRWCellTalus::dGetTotCoarseAdded(void) const
{
   return m_dTotCoarseInput;
}

//! Returns the value for total coarse talus lost
double CRWCellTalus::dGetTotCoarseLost(void) const
{
   return m_dTotCoarseLost;
}
