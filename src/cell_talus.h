/*!
   \class CRWCellTalus
   \brief Real-world class used to represent the talus (unconsolidated sediment resulting from cliff collapse) on a cell layer object
   \details TODO 001 This is a more detailed description of the CRWCellTalus class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cell_sediment.h
   \brief Contains CRWCellTalus definitions
*/

#ifndef TALUS_H
#define TALUS_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
class CRWCellTalus
{
 private:
   //! Current depth equivalent of talus sand sediment in m
   double m_dSand;

   //! Current depth equivalent of talus coarse sediment in m
   double m_dCoarse;

   //! Depth equivalent (m) of talus sand sediment lost this iteration
   double m_dSandLostThisIter;

   //! Depth equivalent (m) of talus coarse sediment lost this iteration
   double m_dCoarseLostThisIter;

   //! Depth equivalent (m) of talus sand sediment added this iteration
   double m_dSandInputThisIter;

   //! Depth equivalent (m) of talus coarse sediment added this iteration
   double m_dCoarseInputThisIter;

   //! Depth equivalent (m) of talus sand sediment added since start of simulation
   double m_dTotSandInput;

   //! Depth equivalent (m) of talus coarse sediment added since start of simulation
   double m_dTotCoarseInput;

   //! Depth equivalent (m) of talus sand sediment lost since start of simulation
   double m_dTotSandLost;

   //! Depth equivalent (m) of talus coarse sediment lost since start of simulation
   double m_dTotCoarseLost;

 protected:
 public:
   CRWCellTalus(void);

   void SetSandDepth(double const);
   double dGetSandDepth(void) const;
   void AddSandDepth(double const);

   void SetCoarseDepth(double const);
   double dGetCoarseDepth(void) const;
   void AddCoarseDepth(double const);

   double dGetSandAddedThisIter(void) const;
   double dGetTotSandAdded(void) const;
   double dGetSandLostThisIter(void) const;
   double dGetTotSandLost(void) const;

   double dGetCoarseAddedThisIter(void) const;
   double dGetTotCoarseAdded(void) const;
   double dGetCoarseLostThisIter(void) const;
   double dGetTotCoarseLost(void) const;
};
#endif // TALUS_H
