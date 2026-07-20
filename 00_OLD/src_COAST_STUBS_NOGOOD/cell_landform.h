/*!
   \class CRWCellLandform
   \brief Real-world class used to represent the landform of a cell
   \details TODO 001 This is a more detailed description of the CCRWCellLandform class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
   \file cell_landform.h
   \brief Contains CRWCellLandform definitions
*/

#ifndef CELL_LANDFORM_H
#define CELL_LANDFORM_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
class CRWCellLandform
{
 private:
   //! Landform category for this cell
   int m_nCategory;

   //! Coast on which this landform sits (if any)
   int m_nCoast;

   //! Point on coast on which this landform sits (if any)
   int m_nPointOnCoastline;

   //! Accumulate wave energy for this landform on this cell
   double m_dAccumWaveEnergy;

   //! The m_uLFData will hold landform data: currently, only cliffs are considered
   union
   {
      struct
      {
         //! Currently unused for all landforms except cliffs
         int m_nDummy;
      } m_sBeachData;

      struct
      {
         //! Cliff notch base elevation
         double m_dNotchApexElev;

         //! Cliff notch, horizontal depth of incision
         double m_dNotchIncision;

         //! Timestep at which cliff collapse occurred
         unsigned long m_ulCollapseTimestep;
      } m_sCliffData;

   } m_uLFData;

 protected:
 public:
   CRWCellLandform();
   ~CRWCellLandform(void);

   void SetLandformCategory(int const);
   int nGetLandformCategory(void) const;

   void SetCoast(int const);
   int nGetCoast(void) const;

   void SetPointOnCoast(int const);
   int nGetPointOnCoast(void) const;

   void SetAccumWaveEnergy(double const);
   void AddToAccumWaveEnergy(double const);
   double dGetAccumWaveEnergy(void) const;

   void SetCliffNotchApexElev(double const);
   double dGetCliffNotchApexElev(void) const;
   void SetCliffNotchIncisionDepth(double const);
   void AddToCliffNotchIncisionDepth(double const);
   double dGetCliffNotchIncisionDepth(void) const;
   // void SetCliffRemaining(double const);
   // double dGetCliffRemaining(void) const;

   void SetCliffCollapseTimestep(unsigned long const);
   unsigned long ulGetCliffCollapseTimestep(void) const;
};
#endif // CELL_LANDFORM_H
