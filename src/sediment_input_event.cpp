/*!
   \file sediment_input_event.cpp
   \brief CRWSedInputEvent routines
   \details TODO 001 A more detailed description of these routines.
   \author Manuel Cobos Budia
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

#include "cme.h"
#include "sediment_input_event.h"

//! Constructor with eight parameters and an initialisation list
CRWSedInputEvent::CRWSedInputEvent(int const nIDIn, unsigned long const ulTimeStepIn, double const dFineIn, double const dSandIn, double const dCoarseIn, double const dLenIn, double const dWidthIn) : //, double const dThickIn):
   m_nLocationID(nIDIn),
   m_ulEventTimeStep(ulTimeStepIn),
   m_dFineSedVol(dFineIn),
   m_dSandSedVol(dSandIn),
   m_dCoarseSedVol(dCoarseIn),
   m_dLen(dLenIn),
   m_dWidth(dWidthIn)
// m_dThick(dThickIn)
{
}

//! Destructor
CRWSedInputEvent::~CRWSedInputEvent(void)
{
}

//! Returns the location ID
int CRWSedInputEvent::nGetLocationID(void) const
{
   return m_nLocationID;
}

//! Returns the timestep of the sediment input event
unsigned long CRWSedInputEvent::ulGetEventTimeStep(void) const
{
   return m_ulEventTimeStep;
}

//! Returns the volume of fine sediment in the sediment input event
double CRWSedInputEvent::dGetFineSedVol(void) const
{
   return m_dFineSedVol;
}

//! Returns the volume of sand sediment in the sediment input event
double CRWSedInputEvent::dGetSandSedVol(void) const
{
   return m_dSandSedVol;
}

//! Returns the volume of coarse sediment in the sediment input event
double CRWSedInputEvent::dGetCoarseSedVol(void) const
{
   return m_dCoarseSedVol;
}

//! Returns the length (in metres?) of the sediment input event
double CRWSedInputEvent::dGetLen(void) const
{
   return m_dLen;
}

//! Returns the width (in metres?) of the sediment input event
double CRWSedInputEvent::dGetWidth(void) const
{
   return m_dWidth;
}

// double CRWSedInputEvent::dGetThick(void)
// {
// return m_dThick;
// }
