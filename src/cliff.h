/*!
   \class CRWCliff
   \brief Real-world class used to represent the 'cliff' category of coastal landform object
   \details TODO 001 This is a more detailed description of the CRWCliff class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cliff.h
   \brief Contains CRWCliff definitions
*/

#ifndef CLIFF_H
#define CLIFF_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "coast.h"
#include "coast_landform.h"

class CRWCliff : public CACoastLandform
{
 private:
   //! Switch to say whether the cliff has just collapsed, earlier in this timestep
   bool m_bCliffHasCollapsed;

   //! The maximum horizontal incision (in external CRS units) of an erosional notch, this is equal to the grid's m_dCellSide
   double m_dMaxNotchIncision;

   //! The horizontal depth (in external CRS units) of the erosional notch, measured inland from the side of the cell that touches the sea
   double m_dNotchIncision;

   //! Elevation (in external CRS units) of the apex of the notch (in CoastalME, the mid point between the top and bottom of the notch). The notch is assumed to extend across the whole width of the coast cell, along the side of the cell that touches the sea
   double m_dNotchApexElev;

 protected:
 public:
   CRWCliff(CRWCoast *, int const, int const, double const, double const, double const, double const);
   ~CRWCliff(void) override;

   void SetCliffCollapsed(void);
   bool bHasCollapsed(void) const;

   void SetNotchApexElev(double const);
   double dGetNotchApexElev(void) const;
   void SetNotchIncision(double const);
   double dGetNotchIncision(void) const;
   void IncreaseNotchIncision(double const);

   bool bReadyToCollapse(double const) const;

   void Display(void) override;
};
#endif // CLIFF_H
