/*!
   \file cell_layer.cpp
   \brief CRWCellLayer routines
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
#include <cstddef>

#include "cell_layer.h"
#include "cell_sediment.h"
#include "cell_talus.h"

//! Constructor
CRWCellLayer::CRWCellLayer(void)
{
   m_pTalus = NULL;
}

//! Destructor
CRWCellLayer::~CRWCellLayer(void)
{
   if (m_pTalus != NULL)
      delete m_pTalus;
}

//! Returns a pointer to the cell's unconsolidated sediment object
CRWCellSediment* CRWCellLayer::pGetUnconsolidatedSediment(void)
{
   return &m_UnconsolidatedSediment;
}

//! Returns a pointer to the cell's consolidated sediment object
CRWCellSediment* CRWCellLayer::pGetConsolidatedSediment(void)
{
   return &m_ConsolidatedSediment;
}

//! Return true if the layer has unconsolidated sediment, false otherwise
bool CRWCellLayer::bHasUncons(void)
{
   if (this->dGetAllUnconsDepth() > 0)
      return true;

   return false;
}

//! Returns the thickness of this cell's fine unconsolidated sediment
double CRWCellLayer::dGetUnconsFineDepth(void) const
{
   return m_UnconsolidatedSediment.dGetFineDepth();
}

//! Returns the thickness of this cell's fine consolidated sediment
double CRWCellLayer::dGetConsFineDepth(void) const
{
   return m_ConsolidatedSediment.dGetFineDepth();
}

//! Returns the thickness of this cell's sand unconsolidated sediment
double CRWCellLayer::dGetUnconsSandDepth(void) const
{
   return m_UnconsolidatedSediment.dGetSandDepth();
}

//! Returns the thickness of this cell's sand consolidated sediment
double CRWCellLayer::dGetConsSandDepth(void) const
{
   return m_ConsolidatedSediment.dGetSandDepth();
}

//! Returns the thickness of this cell's coarse unconsolidated sediment
double CRWCellLayer::dGetUnconsCoarseDepth(void) const
{
   return m_UnconsolidatedSediment.dGetCoarseDepth();
}

//! Returns the thickness of this cell's coarse consolidated sediment
double CRWCellLayer::dGetConsCoarseDepth(void) const
{
   return m_ConsolidatedSediment.dGetCoarseDepth();
}

//! Returns the thickness of this cell's unconsolidated sediment (total for all size classes)
double CRWCellLayer::dGetAllUnconsDepth(void) const
{
   return (m_UnconsolidatedSediment.dGetFineDepth() + m_UnconsolidatedSediment.dGetSandDepth() + m_UnconsolidatedSediment.dGetCoarseDepth());
}

//! Returns the thickness of this cell's consolidated sediment (total for all size classes)
double CRWCellLayer::dGetAllConsDepth(void) const
{
   return (m_ConsolidatedSediment.dGetFineDepth() + m_ConsolidatedSediment.dGetSandDepth() + m_ConsolidatedSediment.dGetCoarseDepth());
}

//! Returns the thickness of this cell's sediment (total for all size classes, both consolidated and unconsolidated)
double CRWCellLayer::dGetTotalThickness(void) const
{
   return (m_UnconsolidatedSediment.dGetFineDepth() + m_UnconsolidatedSediment.dGetSandDepth() + m_UnconsolidatedSediment.dGetCoarseDepth() + m_ConsolidatedSediment.dGetFineDepth() + m_ConsolidatedSediment.dGetSandDepth() + m_ConsolidatedSediment.dGetCoarseDepth());
}

// //! Returns the thickness of unconsolidated sediment lost due to notch incision (total for all size classes)
// double CRWCellLayer::dGetNotchUnconsolidatedLost(void) const
// {
// return (m_UnconsolidatedSediment.dGetNotchFineLost() + m_UnconsolidatedSediment.dGetNotchSandLost() + m_UnconsolidatedSediment.dGetNotchCoarseLost());
// }

// //! Returns the thickness of consolidated sediment lost due to notch incision (total for all size classes)
// double CRWCellLayer::dGetNotchConsolidatedLost(void) const
// {
// return (m_ConsolidatedSediment.dGetNotchFineLost() + m_ConsolidatedSediment.dGetNotchSandLost() + m_ConsolidatedSediment.dGetNotchCoarseLost());
// }

// double CRWCellLayer::dGetVolSedFraction(void) const
// {
// return m_VdolSedFraction;
// }

// void CRWCellLayer::SetVolSedFraction(double const dNewVolSedFraction)
// {
// m_VdolSedFraction = dNewVolSedFraction;
// }
//
// double CRWCellLayer::dGetMechResistance(void) const
// {
// return m_dMechResistance;
// }

// void CRWCellLayer::SetMechResistance(double const dNewMechResistance)
// {
// m_dMechResistance = dNewMechResistance;
// }

// double CRWCellLayer::dGetConsolidationStatus(void) const
// {
// return m_dConsolidationStatus;
// }

// void CRWCellLayer::SetConsolidationStatus(double const dNewConsolidationStatus)
// {
// m_dConsolidationStatus = dNewConsolidationStatus;
// }

//! Returns true if the layer has talus, false otherwise
bool CRWCellLayer::bHasTalus(void)
{
   if (m_pTalus == NULL)
      return false;

   return true;
}

//! Returns a pointer to the layer's talus object, or NULL if there is no talus object
CRWCellTalus* CRWCellLayer::pGetTalus(void) const
{
   return m_pTalus;
}

//! Returns a pointer to the layer's talus object. If there is no talus object, then create one
CRWCellTalus* CRWCellLayer::pGetOrCreateTalus(void)
{
   if (m_pTalus == NULL)
   {
      // No talus object here, so create one
      m_pTalus = new CRWCellTalus();
   }

   return m_pTalus;
}

// Deletes this layer's talus object and resets the pointer
void CRWCellLayer::DeleteTalus(void)
{
   delete m_pTalus;
   m_pTalus = NULL;
}

