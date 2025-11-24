/*!
   \class CRWCellSediment
   \brief Real-world class used to represent the sediment (either consolidated or unconsolidated) associated with a cell layer object
   \details TODO 001 This is a more detailed description of the CRWCellSediment class.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
   \file cell_sediment.h
   \brief Contains CRWCellSediment definitions
*/

#ifndef SEDIMENT_H
#define SEDIMENT_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
class CRWCellSediment
{
 private:
   //! Depth equivalent of fine sediment in m
   double m_dFine;

   //! Depth equivalent (m) of fine sediment lost via notch incision
   double m_dNotchFineLost;

   //! Depth equivalent of sand sediment in m
   double m_dSand;

   //! Depth equivalent (m) of sand sediment lost via notch incision
   double m_dNotchSandLost;

   //! Depth equivalent of coarse sediment in m
   double m_dCoarse;

   //! Depth equivalent (m) of coarse sediment lost via notch incision
   double m_dNotchCoarseLost;

   //! Depth equivalent (m) of fine sediment added via sediment input events, this iteration
   double m_dFineSedimentInputThisIter;

   //! Depth equivalent (m) of sand sediment added via sediment input events, this iteration
   double m_dSandSedimentInputThisIter;

   //! Depth equivalent (m) of coarse sediment added via sediment input events, this iteration
   double m_dCoarseSedimentInputThisIter;

   //! Depth equivalent (m) of fine sediment added via sediment input events, since start of simulation
   double m_dTotFineSedimentInput;

   //! Depth equivalent (m) of sand sediment added via sediment input events, since start of simulation
   double m_dTotSandSedimentInput;

   //! Depth equivalent (m) of coarse sediment added via sediment input events, since start of simulation
   double m_dTotCoarseSedimentInput;

 protected:
 public:
   CRWCellSediment(void);
   CRWCellSediment(CRWCellSediment const&); // Copy constructor defined explicitly, to stop cppcheck from complaining

   CRWCellSediment& operator=(const CRWCellSediment&);

   void SetFineDepth(double const);
   double dGetFineDepth(void) const;
   void AddFineDepth(double const);

   void SetSandDepth(double const);
   double dGetSandDepth(void) const;
   void AddSandDepth(double const);

   void SetCoarseDepth(double const);
   double dGetCoarseDepth(void) const;
   void AddCoarseDepth(double const);

   void SetNotchFineLost(double const);
   // void IncrNotchFineLost(double const);
   double dGetNotchFineLost(void) const;

   void SetNotchSandLost(double const);
   // void IncrNotchSandLost(double const);
   double dGetNotchSandLost(void) const;

   void SetNotchCoarseLost(double const);
   // void IncrNotchCoarseLost(double const);
   double dGetNotchCoarseLost(void) const;

   void AddFineSedimentInputDepth(double const);
   void AddSandSedimentInputDepth(double const);
   void AddCoarseSedimentInputDepth(double const);
   double dGetFineSedimentInputDepth(void) const;
   double dGetSandSedimentInputDepth(void) const;
   double dGetCoarseSedimentInputDepth(void) const;
   double dGetTotAllSedimentInputDepth(void) const;
   void InitThisIterSedimentInputAll(void);
};
#endif // SEDIMENT_H
