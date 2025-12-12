/*!
   \class CRWSedInputEvent
   \brief Class used to represent a sediment input event
   \details This class represent a sediment input event such as sediment derived from inland (e.g. at the mouth of a gully or rambla), or sediment from an intervention such as beach nourishment
   \author Manuel Cobos Budia
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file sediment_input_event.h
   \brief Contains CRWSedInputEvent definitions
*/

#ifndef CSEDINPUT_H
#define CSEDINPUT_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
class CRWSedInputEvent
{
 private:
   //! The location ID in the shapefile
   int m_nLocationID;

   //! The timing of the sediment input event
   unsigned long m_ulEventTimeStep;

   //! The volume (m3) of fine sediment in the sediment input event
   double m_dFineSedVol;

   //! The volume (m3) of sand sediment in the sediment input event
   double m_dSandSedVol;

   //! The volume (m3) of coarse sediment in the sediment input event
   double m_dCoarseSedVol;

   //! The coast-normal length (m) of the sediment block
   double m_dLen;

   //! The along-coast width (m) of the sediment block
   double m_dWidth;

 public:
   CRWSedInputEvent(int const, unsigned long const, double const, double const, double const, double const, double const); //, double const);
   ~CRWSedInputEvent(void);

   int nGetLocationID(void) const;
   unsigned long ulGetEventTimeStep(void) const;
   double dGetFineSedVol(void) const;
   double dGetSandSedVol(void) const;
   double dGetCoarseSedVol(void) const;
   double dGetLen(void) const;
   double dGetWidth(void) const;
   // double dGetThick(void);
};
#endif // CSEDINPUT_H
