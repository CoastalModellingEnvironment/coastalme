/*!
   \file cell_landform.cpp
   \brief CRWCellLandform routines
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
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
#include "cell_landform.h"

//! Constructor
CRWCellLandform::CRWCellLandform()
    : m_nCategory(LF_UNKNOWN),
      m_nCoast(-1),
      m_nPointOnCoastline(-1),
      m_dAccumWaveEnergy(0)
{
   m_uLFData.m_sCliffData.m_dNotchApexElev = DBL_NODATA;
   m_uLFData.m_sCliffData.m_dNotchIncision = DBL_NODATA;
   m_uLFData.m_sCliffData.m_ulCollapseTimestep = UNSIGNED_LONG_NODATA;
}

//! Destructor
CRWCellLandform::~CRWCellLandform(void)
{
}

//! Set the landform category
void CRWCellLandform::SetLFCategory(int const nClassIn)
{
   m_nCategory = nClassIn;
}

//! Get the landform category
int CRWCellLandform::nGetLFCategory(void) const
{
   return m_nCategory;
}

//! Set the coast number
void CRWCellLandform::SetCoast(int const nCoastIn)
{
   m_nCoast = nCoastIn;
}

//! Get the coast number
int CRWCellLandform::nGetCoast(void) const
{
   return m_nCoast;
}

//! Set the number of the point on the coastline
void CRWCellLandform::SetPointOnCoast(int const nPointOnCoastIn)
{
   m_nPointOnCoastline = nPointOnCoastIn;
}

//! Set the number of the point on the coastline
int CRWCellLandform::nGetPointOnCoast(void) const
{
   return m_nPointOnCoastline;
}

//! Set accumulated wave energy
void CRWCellLandform::SetAccumWaveEnergy(double const dEnergyIn)
{
   m_dAccumWaveEnergy = dEnergyIn;
}

//! Add to acccumulated wave energy
void CRWCellLandform::AddToAccumWaveEnergy(double const dEnergyIn)
{
   m_dAccumWaveEnergy += dEnergyIn;
}

//! Get accumulated wave energy
double CRWCellLandform::dGetAccumWaveEnergy(void) const
{
   return m_dAccumWaveEnergy;
}

//! Set cliff notch apex elevation
void CRWCellLandform::SetCliffNotchApexElev(double const dElevIn)
{
   m_uLFData.m_sCliffData.m_dNotchApexElev = dElevIn;
}

//! Get cliff notch apex elevation: this is the elevation of the most deeply incised horizontal plane within the notch, see Figure 3 in Trenhaile, A.S. (2015). Coastal notches: Their morphology, formation, and function. Earth-Science Reviews 150, 285-304. In CoastalME, where the notch is assumed to have the same depth of incision on all horizonal planes within the notch, the apex is assumed to be the vertical mid-point between the upper and lower edges of the notch
double CRWCellLandform::dGetCliffNotchApexElev(void) const
{
   return m_uLFData.m_sCliffData.m_dNotchApexElev;
}

//! Set the horizontal depth of cliff notch incision for this cell. This is the horizontal distance from the edge of the cell to the incised back of the notch
void CRWCellLandform::SetCliffNotchIncisionDepth(double const dLenIn)
{
   m_uLFData.m_sCliffData.m_dNotchIncision = dLenIn;
}

//! Adds to the horizontal depth of cliff notch incision for this cell. This is the horizontal distance from the edge of the cell to the incised back of the notch
void CRWCellLandform::AddToCliffNotchIncisionDepth(double const dLenIn)
{
   m_uLFData.m_sCliffData.m_dNotchIncision += dLenIn;
}

//! Get the horizontal depth of cliff notch incision for this cell. This is the horizontal distance from the edge of the cell to the incised back of the notch
double CRWCellLandform::dGetCliffNotchIncisionDepth(void) const
{
   return m_uLFData.m_sCliffData.m_dNotchIncision;
}

// //! Set the cliff depth remaining on this cell
// void CRWCellLandform::SetCliffRemaining(double const dLenIn)
// {
//    m_uLFData.m_sCliffData.m_dRemaining = dLenIn;
// }

//! Set the timestep at which cliff collapse occurred
void CRWCellLandform::SetCliffCollapseTimestep(unsigned long const ulTimestep)
{
   m_uLFData.m_sCliffData.m_ulCollapseTimestep = ulTimestep;
}

//! Returns the timestep at which cliff collapse occurred
unsigned long CRWCellLandform::ulGetCliffCollapseTimestep(void) const
{
   return m_uLFData.m_sCliffData.m_ulCollapseTimestep;
}

