/*!
   \class CRWCellLayer
   \brief Real-world class used to represent the sediment layers associated with a cell object
   \details TODO 001 This is a more detailed description of the CCRWCellLayer class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cell_layer.h
   \brief Contains CRWCellLayer definitions
*/

#ifndef CELL_LAYER_H
#define CELL_LAYER_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "cell_sediment.h"
#include "cell_talus.h"

class CRWCellLayer
{
 private:
   // double m_VdolSedFraction;
   // double m_dMechResistance;
   // double m_dConsolidationStatus;

   //! This cell's unconsolidated sediment object
   CRWCellSediment m_UnconsolidatedSediment;

   //! This cell's consolidated sediment object
   CRWCellSediment m_ConsolidatedSediment;

   //! Pointer to this cell layer's talus (unconsolidated sediment derived from cliff collapse), if any
   CRWCellTalus* m_pTalus;

 protected:
 public:
   CRWCellLayer(void);
   ~CRWCellLayer(void);

   CRWCellSediment* pGetUnconsolidatedSediment(void);
   CRWCellSediment* pGetConsolidatedSediment(void);

   bool bHasUncons(void);
   double dGetUnconsFineDepth(void) const;
   double dGetConsFineDepth(void) const;
   double dGetUnconsSandDepth(void) const;
   double dGetConsSandDepth(void) const;
   double dGetUnconsCoarseDepth(void) const;
   double dGetConsCoarseDepth(void) const;

   double dGetAllUnconsDepth(void) const;
   double dGetAllConsDepth(void) const;
   double dGetTotalThickness(void) const;

   // double dGetNotchUnconsolidatedLost(void) const;
   // double dGetNotchConsolidatedLost(void) const;

   // double dGetVolSedFraction(void) const;
   void SetVolSedFraction(double const);
   // double dGetMechResistance(void) const;
   // void SetMechResistance(double const);
   // double dGetConsolidationStatus(void) const;
   // void SetConsolidationStatus(double const);

   bool bHasTalus(void);
   CRWCellTalus* pGetTalus(void) const;
   CRWCellTalus* pGetOrCreateTalus(void);
   void DeleteTalus(void);
};
#endif // CELL_LAYER_H
