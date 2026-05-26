/*!
   \file do_cliff_collapse.cpp
   \brief Collapses cliffs if a critical notch depth is exceeded
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <cstddef>
#include <assert.h>

#include <iostream>
// using std::cerr;
using std::endl;
using std::ios;

#include "cme.h"
#include "simulation.h"
#include "cliff.h"
#include "cell_talus.h"
#include "coast_landform.h"
#include "2di_point.h"

//===============================================================================================================================
//! Update accumulated wave energy in coastal landform objects. If the object is a cliff, then deepen the incised notch. If the notch is sufficiently deep, cliff collapse occurs.
//! CoastalME's representation of notch incision is based on Trenhaile, A.S. (2015). Coastal notches: Their morphology, formation, and function. Earth-Science Reviews 150, 285-304
//===============================================================================================================================
int CSimulation::nDoAllWaveEnergyToCoastLandforms(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Calculating cliff collapse" << endl;

   int nRet;

   // First go along each coastline and at each point on the coastline, update the total wave energy which it has experienced TODO Note that currently, only cliff objects respond to accumulated wave energy
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint);

         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

         // TODO Talus is protecting this cell
         if (m_pRasterGrid->m_Cell[nX][nY].dGetTalusDepth() > 0)
            continue;

         // First get wave energy for the coastal landform object
         double const dWaveHeightAtCoast = m_VCoast[nCoast].dGetCoastWaveHeight(nCoastPoint);

         // If the waves at this point are off-shore, then do nothing, just move to next coast point
         if (bFPIsEqual(dWaveHeightAtCoast, DBL_NODATA, TOLERANCE))
            continue;

         // OK we have on-shore waves so get the previously-calculated wave energy
         double const dWaveEnergy = m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nCoastPoint);

         // And save the accumulated value
         pCoastLandform->IncTotAccumWaveEnergy(dWaveEnergy);

         int const nCat = pCoastLandform->nGetLandFormCategory();

         // Is this a cliff?
         if (nCat == LF_CLIFF)
         {
            // It is, so get the cliff object
            CRWCliff* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

            // And do the notch incision, if any. Note that we consider sediment eroded due to notch incision to be still in place until cliff collapse, i.e. the sediment which filled the notch, pre-incision, is assumed to remain there. If the notch is eventually incised sufficiently to cause cliff collapse, then the sediment from the notch volume is included with the above-notch talus
            if (! bIncreaseCliffNotchIncision(nCoast, nX, nY, pCliff, dWaveEnergy))
               // No incision of this cliff
               continue;

            // OK, we've had some incision. So is the notch now extended enough to cause collapse (either because the overhang is greater than the threshold overhang, or because there is no sediment remaining)?
            if (pCliff->bReadyToCollapse(m_dNotchIncisionAtCollapse))
            {
               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t ready to collapse, notch incision = " << pCliff->dGetNotchIncision() << " threshold incision = " << m_dNotchIncisionAtCollapse << endl;

               // // DEBUG CODE ============================================================================================================================================
               // // Get total depths of sand consolidated and unconsolidated for every cell
               // if (m_ulIter == 5)
               // {
               // double dTmpSandCons = 0;
               // double dTmpSandUncons = 0;
               // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
               // {
               // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
               // {
               // dTmpSandCons += m_pRasterGrid->m_Cell[nX1][nY1].dGetConsSandDepthAllLayers();
               //
               // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
               // }
               // }
               //
               //    // Get the cliff cell's grid coords
               // int nXCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
               // int nYCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();
               //
               //    // Get this cell's polygon
               // int nPoly = m_pRasterGrid->m_Cell[nXCliff][nYCliff].nGetPolygonID();
               //
               // LogStream << endl;
               // LogStream << "*****************************" << endl;
               // LogStream << m_ulIter << ": before cliff collapse on nPoly = " << nPoly << " total consolidated sand = " << dTmpSandCons * m_dCellArea << " total unconsolidated sand = " << dTmpSandUncons * m_dCellArea << endl;
               // }
               // // DEBUG CODE ============================================================================================================================================

               // It is ready to collapse
               int nNotchLayer;
               double dCliffElevPreCollapse = 0;
               double dCliffElevPostCollapse = 0;
               double dFineCollapse = 0;
               double dSandCollapse = 0;
               double dCoarseCollapse = 0;

               // Do the cliff collapse
               nRet = nDoCliffCollapse(nCoast, pCliff, dFineCollapse, dSandCollapse, dCoarseCollapse, nNotchLayer, dCliffElevPreCollapse, dCliffElevPostCollapse);
               if (nRet != RTN_OK)
               {
                  if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     LogStream << m_ulIter << ":\t" << WARN << "problem with cliff collapse, continuing however" << endl;
               }

               // Deposit all sand and/or coarse sediment derived from this cliff collapse as talus, on the cell on which collapse occurred
               DoCliffCollapseTalusDeposition(/*nCoast,*/ pCliff, dSandCollapse, dCoarseCollapse, nNotchLayer);

               // // DEBUG CODE ============================================================================================================================================
               // // Get total depths of sand consolidated and unconsolidated for every cell
               // if (m_ulIter == 5)
               // {
               // double dTmpSandCons = 0;
               // double dTmpSandUncons = 0;
               // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
               // {
               // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
               // {
               // dTmpSandCons += m_pRasterGrid->m_Cell[nX1][nY1].dGetConsSandDepthAllLayers();
               //
               // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
               // }
               // }
               //
               //    // Get the cliff cell's grid coords
               // int nXCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
               // int nYCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();
               //
               //    // Get this cell's polygon
               // int nPoly = m_pRasterGrid->m_Cell[nXCliff][nYCliff].nGetPolygonID();
               //
               // LogStream << endl;
               // LogStream << "*****************************" << endl;
               // LogStream << m_ulIter << ": after cliff collapse on nPoly = " << nPoly << " total consolidated sand = " << dTmpSandCons * m_dCellArea << " total unconsolidated sand = " << dTmpSandUncons * m_dCellArea << endl;
               // LogStream << m_ulIter << ": total consolidated sand lost this iteration =  " << (m_dStartIterConsSandAllCells - dTmpSandCons) * m_dCellArea << endl;
               // LogStream << m_ulIter << ": total unconsolidated sand added this iteration =  " << (dTmpSandUncons - m_dStartIterUnconsSandAllCells) * m_dCellArea << endl;
               //
               // double dTmpAllPolySandErosion = 0;
               // double dTmpAllPolySandDeposition = 0;
               // for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
               // {
               // double dTmpSandErosion = m_pVCoastPolygon[n]->dGetCliffCollapseErosionSand() * m_dCellArea ;
               // double dTmpSandDeposition = m_pVCoastPolygon[n]->dGetCliffCollapseUnconsSandDeposition() * m_dCellArea ;
               //
               // LogStream << m_ulIter << ": polygon = " << m_pVCoastPolygon[n]->nGetPolygonThisCoastID() << " sand erosion = " << dTmpSandErosion << " sand deposition = " << dTmpSandDeposition << endl;
               //
               // dTmpAllPolySandErosion += dTmpSandErosion;
               // dTmpAllPolySandDeposition += dTmpSandDeposition;
               // }
               //
               // LogStream << "-------------------------------------------" << endl;
               // LogStream << m_ulIter << ": all polygons, sand erosion = " << dTmpAllPolySandErosion << " sand deposition = " << dTmpAllPolySandDeposition << endl;
               // LogStream << "*****************************" << endl;
               // }
               // // DEBUG CODE ============================================================================================================================================
            }
         }
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_ALL)
      LogStream << m_ulIter << ":\t total cliff collapse (m^3) = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons + m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << " (fine = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons) * m_dCellArea << ", sand = " << (m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons) * m_dCellArea << ", coarse = " << (m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << "), talus deposition (m^3) = " << (m_dThisIterUnconsSandCliffDeposition + m_dThisIterUnconsCoarseCliffDeposition) * m_dCellArea << " (sand = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ", coarse = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ")" << endl << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Simulates cliff collapse on a single cell. Collapse happens when when a notch which is incised into the cell's consolidated sediment layer exceeds a critical horizontal incision. This routine updates the cliff object, the cell 'under' the cliff object, and the polygon which contains the cliff object
//===============================================================================================================================
int CSimulation::nDoCliffCollapse(int const nCoast, CRWCliff* pCliff, double& dFineCollapse, double& dSandCollapse, double& dCoarseCollapse, int& nNotchLayer, double& dPreCollapseCellElev, double& dPostCollapseCellElevNoTalus)
{
   // Get the cliff cell's grid coords
   int const nX = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
   int const nY = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();

   // Get this cell's polygon
   int const nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
   if (nPoly == INT_NODATA)
   {
      // This cell isn't in a polygon
      LogStream << m_ulIter << ":\t" << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} is not in a polygon" << endl;
      return RTN_ERR_CLIFF_NOT_IN_POLYGON;
   }

   // Get a pointer to the polygon
   CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   // Get the elevation of the apex of the notch from the cliff object
   double const dNotchElev = pCliff->dGetNotchApexElev();

   // Get the index of the layer containing the notch (layer 0 being just above basement)
   nNotchLayer = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dNotchElev);
   double const dTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevOmitTalus();
   double const dTopElevIncTalus = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevIncTalus();

   // Safety check: is the notch elevation above the top of the consolidated sediment? If so, do no more
   if (nNotchLayer == ELEV_ABOVE_SEDIMENT_TOP)
   {
      LogStream << m_ulIter << ": cliff ready to collapse at [" << nX << "][" << nY << "]  = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} nNotchLayer is above sediment top, notch layer = " << nNotchLayer << " dNotchElev = " << dNotchElev << " m_dNotchApexAboveMHW = " << m_dNotchApexAboveMHW << " dPreCollapseCellElev = " << dPreCollapseCellElev << " cons sediment top elev without talus = " << dTopElevNoTalus << " cons sediment top elev inc talus = " << dTopElevIncTalus << endl;

      return RTN_OK;
   }

   // More safety checks
   if (nNotchLayer == ELEV_IN_BASEMENT)
   {
      LogStream << m_ulIter << ":\t" << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} nNotchLayer is in basement, notch layer = " << nNotchLayer << " dNotchElev = " << dNotchElev << " m_dNotchApexAboveMHW = " << m_dNotchApexAboveMHW << " dPreCollapseCellElev = " << dPreCollapseCellElev << " cons sediment top elev without talus = " << dTopElevNoTalus << " cons sediment top elev inc talus = " << dTopElevIncTalus << endl;
      return RTN_ERR_CLIFF_NOTCH;
   }

   if (nNotchLayer < 0)
   {
      LogStream << m_ulIter << ":\t" << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} nNotchLayer less than zero, notch layer = " << nNotchLayer << " dNotchElev = " << dNotchElev << " m_dNotchApexAboveMHW = " << m_dNotchApexAboveMHW << " dPreCollapseCellElev = " << dPreCollapseCellElev << " cons sediment top elev without talus = " << dTopElevNoTalus << " cons sediment top elev inc talus = " << dTopElevIncTalus << endl;
      return RTN_ERR_CLIFF_NOTCH;
   }

   // Notch layer is OK, so flag the coastline cliff object as having collapsed
   pCliff->SetCliffCollapsed();

   int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();

   // Safety check
   if (nTopLayer == INT_NODATA)
   {
      LogStream << m_ulIter << ":\t" << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] nTopLayer = " << nTopLayer << endl;
      return RTN_ERR_NO_TOP_LAYER_DURING_CLIFF_COLLAPSE_CALC;
   }

   // Set flags to say that the notch layer, and all layers above it, have changed
   for (int nLayer = nNotchLayer; nLayer <= nTopLayer; nLayer++)
   {
      m_bConsSedChangedThisIter[nLayer] = true;
      m_bUnconsChangedThisIter[nLayer] = true;
   }

   // Get the pre-collapse cliff elevation (we assume that there is no talus on top of the cliff)
   dPreCollapseCellElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

   // Now calculate the vertical depth of sediment lost in this cliff collapse, note that this includes the sediment which filled the notch before any incision took place
   double dAvailable = 0;
   double dFineConsLost = 0;
   double dFineUnconsLost = 0;
   double dSandConsLost = 0;
   double dSandUnconsLost = 0;
   double dCoarseConsLost = 0;
   double dCoarseUnconsLost = 0;

   // Update the cell's sediment. If there are sediment layers above the notched layer, we must remove sediment from the whole depth of each layer
   for (int n = nTopLayer; n > nNotchLayer; n--)
   {
      // Start with the unconsolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetFineDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchFineLost();

      if (dAvailable > 0)
      {
         dFineCollapse += dAvailable;
         dFineUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetFineDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchFineLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetSandDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchSandLost();

      if (dAvailable > 0)
      {
         dSandCollapse += dAvailable;
         dSandUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetSandDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchSandLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetCoarseDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();

      if (dAvailable > 0)
      {
         dCoarseCollapse += dAvailable;
         dCoarseUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetCoarseDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchCoarseLost(0);
      }

      // Now get the consolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetFineDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchFineLost();

      if (dAvailable > 0)
      {
         dFineCollapse += dAvailable;
         dFineConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetFineDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchFineLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetSandDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchSandLost();

      if (dAvailable > 0)
      {
         dSandCollapse += dAvailable;
         dSandConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetSandDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchSandLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetCoarseDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchCoarseLost();

      if (dAvailable > 0)
      {
         dCoarseCollapse += dAvailable;
         dCoarseConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetCoarseDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchCoarseLost(0);
      }
   }

   // Now calculate the sediment lost from the consolidated layer into which the erosional notch was incised
   double const dNotchLayerTop = m_pRasterGrid->m_Cell[nX][nY].dCalcLayerElev(nNotchLayer);
   double const dNotchLayerThickness = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->dGetTotalThickness();
   double const dNotchLayerFracRemoved = (dNotchLayerTop - dNotchElev) / dNotchLayerThickness;

   // Safety check
   double dTmp = dNotchElev + dFineConsLost + dFineUnconsLost + dSandConsLost + dSandUnconsLost + dCoarseConsLost + dCoarseUnconsLost;
   if (dTmp > m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus())
      LogStream << m_ulIter << ":\t TOO MUCH SEDIMENT AT CLIFF COLLAPSE sediment depth = " << dTmp << " sediment top elevation inc talus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << endl;

   if (dTmp > m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus())
      LogStream << m_ulIter << ":\t TOO MUCH SEDIMENT AT CLIFF COLLAPSE sediment depth = " << dTmp << " sediment top elevation noy inc talus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << endl;

   // Sort out the notched layer's sediment, both consolidated and unconsolidated, for this cell. First the unconsolidated sediment
   double dFineDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
   dAvailable = dFineDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchFineLost();

   if (dAvailable > 0)
   {
      // Some unconsolidated fine sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dFineDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchFineLost(0);
   }

   double dSandDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
   dAvailable = dSandDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchSandLost();

   if (dAvailable > 0)
   {
      // Some unconsolidated sand sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchSandLost(0);
   }

   double dCoarseDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
   dAvailable = dCoarseDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();

   if (dAvailable > 0)
   {
      // Some unconsolidated coarse sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchCoarseLost(0);
   }

   // Do the same for the consolidated sediment
   dFineDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetFineDepth();
   dAvailable = dFineDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchFineLost();

   if (dAvailable > 0)
   {
      // Some consolidated fine sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetFineDepth(dFineDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchFineLost(0);
   }

   dSandDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetSandDepth();
   dAvailable = dSandDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchSandLost();

   if (dAvailable > 0)
   {
      // Some consolidated sand sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetSandDepth(dSandDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchSandLost(0);
   }

   dCoarseDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();
   dAvailable = dCoarseDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchCoarseLost();

   if (dAvailable > 0)
   {
      // Some consolidated coarse sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetCoarseDepth(dCoarseDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchCoarseLost(0);
   }

   // Update the cell's totals for cliff collapse erosion
   m_pRasterGrid->m_Cell[nX][nY].IncrCliffCollapseErosion(dFineCollapse, dSandCollapse, dCoarseCollapse);

   // Update the cell's layer elevations (pre talus deposition) and d50
   m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

   // Get the post-collapse cell elevation (talus will be deposited above this)
   dPostCollapseCellElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      LogStream << m_ulIter << ":\t coast " << nCoast << " cliff collapse at [" << nX << "][" << nY << "]" << ((dSandCollapse > 0) || (dCoarseCollapse > 0) ? ", before talus deposition:" : "") << " original cell elevation = " << dPreCollapseCellElev << ", new cell elevation = " << dPostCollapseCellElevNoTalus << ", change in elevation = " << dPreCollapseCellElev - dPostCollapseCellElevNoTalus << endl;

   // Update this-polygon totals: add to the depths of cliff collapse erosion for this polygon
   pPolygon->AddCliffCollapseErosionFine(dFineCollapse);
   pPolygon->AddCliffCollapseToSuspensionFine(dFineCollapse);
   pPolygon->AddCliffCollapseErosionSand(dSandCollapse);
   pPolygon->AddCliffCollapseErosionCoarse(dCoarseCollapse);

   // And update the this-timestep totals and the grand totals for the number of cells with cliff collapse
   m_nNumThisIterCliffCollapse++;
   m_nNumTotCliffCollapse++;

   // Add to this-iteration totals of fine sediment (consolidated and unconsolidated) eroded via cliff collapse
   m_dThisIterCliffCollapseErosionFineUncons += dFineUnconsLost;
   m_dThisIterCliffCollapseErosionFineCons += dFineConsLost;

   // Also add to the total suspended load. Note that this addition to the suspended load has not yet been added to all cells, this happens in nEndOfTimestepUpdateGrid()
   m_dThisIterFineSedimentToSuspension += (dFineConsLost + dFineUnconsLost);

   // Add to this-iteration totals of sand and coarse sediment (consolidated and unconsolidated) eroded via cliff collapse
   m_dThisIterCliffCollapseErosionSandUncons += dSandUnconsLost;
   m_dThisIterCliffCollapseErosionSandCons += dSandConsLost;
   m_dThisIterCliffCollapseErosionCoarseUncons += dCoarseUnconsLost;
   m_dThisIterCliffCollapseErosionCoarseCons += dCoarseConsLost;

   // Save the timestep at which cliff collapse occurred
   m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffCollapseTimestep(m_ulIter);

   // Reset cell cliff info
   m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchIncisionDepth(m_dCellSide);

   // Final safety check
   int const nNewTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();
   if (nNewTopLayer == INT_NODATA)
      return RTN_ERR_NO_TOP_LAYER_DURING_CLIFF_COLLAPSE_CALC;

   return RTN_OK;
}

//===============================================================================================================================
//! Increase the incision (if any) of a cliff notch, assuming a linear decrease in incision with distance downwards from notch apex. Returns false if no incision
//===============================================================================================================================
bool CSimulation::bIncreaseCliffNotchIncision(int const nCoast, int const nX, int const nY, CRWCliff* pCliff, double const dWaveEnergy)
{
   // Get the coastline point of the cliff
   int const nCoastPoint = pCliff->nGetPointOnCoast();

   // And get the wave runup at this point
   double const dRunup = m_VCoast[nCoast].dGetRunUp(nCoastPoint);
   double const dWaveElev = m_dThisIterSWL + dRunup;

   // Get the apex elevation of the cliff notch
   double const dNotchApexElev = pCliff->dGetNotchApexElev();

   // Is there a notch?
   if (! bFPIsEqual(dNotchApexElev, DBL_NODATA, TOLERANCE))
   {
      // This is a notch in this cliff object, so get the cutoff elevation (if this-iteration SWL is below this, there is no incision)
      double const dCutoffElev = dNotchApexElev - CLIFF_NOTCH_CUTOFF_DISTANCE;

      if (dWaveElev < dCutoffElev)
      {
         // SWL is below the cutoff elevation, so no incision of this existing notch
         LogStream << m_ulIter << ":\t NO incision of existing notch at [" << nX << "][" << nY << "] dWaveElev = " << dWaveElev << " dCutoffElev = " << dCutoffElev << " dRunup = " << dRunup << " dNotchApexElev = " << dNotchApexElev << " Sediment top without talus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus() << " sediment top with talus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << endl;

         return false;
      }

      // SWL is above the cutoff, so we have some incision of this existing notch
      double dWeight;
      if (dWaveElev > dNotchApexElev)
         dWeight = 1;
      else
         // Assume a linear decrease in incision with distance downwards from notch apex
         dWeight = 1 - ((dNotchApexElev - dWaveElev) / CLIFF_NOTCH_CUTOFF_DISTANCE);

      // Calculate this-timestep cliff notch erosion (is a length in external CRS units). Note that only consolidated sediment can have a cliff notch
      double const dNotchIncision = dWeight * dWaveEnergy / m_dCliffErosionResistance;

      // Deepen the cliff object's erosional notch as a result of wave energy during this timestep. Note that notch deepening may be constrained, since this-timestep notch extension cannot exceed the length (i.e. cellside minus notch depth) of sediment remaining on the cell
      pCliff->IncreaseNotchIncision(dNotchIncision);

      // And add to the cell's accumulated wave energy
      m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->AddToAccumWaveEnergy(dWaveEnergy * dWeight);

      LogStream << m_ulIter << ":\t incision of existing notch at [" << nX << "][" << nY << "] accumulated wave energy = " << m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy() << " dWaveElev = " << dWaveElev << " dCutoffElev = " << dCutoffElev << " dRunup = " << dRunup << "  dWeight = " << dWeight << " dNotchApexElev = " << dNotchApexElev << " increase in notch incision = " << dNotchIncision << " total notch incision = " << pCliff->dGetNotchIncision() << " threshold incision = " << m_dNotchIncisionAtCollapse << endl;

      return true;
   }
   else
   {
      // No notch in this cliff object. Can we create one?
      double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

      if (m_dThisIterNewNotchApexElev < dSedTopElevNoTalus)
      {
         // Yes we can create a notch here
         pCliff->SetNotchApexElev(m_dThisIterNewNotchApexElev - SED_ELEV_TOLERANCE);
         pCliff->SetNotchIncision(0);

         // Get the cutoff elevation (if this-iteration SWL is below this, there is no incision)
         double const dCutoffElev = m_dThisIterNewNotchApexElev - CLIFF_NOTCH_CUTOFF_DISTANCE;

         if (dWaveElev < dCutoffElev)
         {
            // SWL is below the cutoff elevation, so no incision of this existing notch
            LogStream << m_ulIter << ":\t NO incision of new notch at [" << nX << "][" << nY << "] dWaveElev = " << dWaveElev << " dCutoffElev = " << dCutoffElev << " dRunup = " << dRunup << " dNotchApexElev = " << dNotchApexElev << " sediment top without talus = " << dSedTopElevNoTalus << " sediment top with talus = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << endl;

            return false;
         }

         // We have some notch incision of this newly-created notch
         double dWeight;
         if (dWaveElev > dNotchApexElev)
            // Should not happen: safety check
            dWeight = 1;
         else
            // Assume a linear decrease in incision with distance downwards from notch apex
            dWeight = 1 - ((dNotchApexElev - dWaveElev) / CLIFF_NOTCH_CUTOFF_DISTANCE);

         // Calculate this-timestep cliff notch erosion (is a length in external CRS units). Note that only consolidated sediment can have a cliff notch
         double const dNotchIncision = dWeight * dWaveEnergy / m_dCliffErosionResistance;

         // Deepen the cliff object's erosional notch as a result of wave energy during this timestep. Note that notch deepening may be constrained, since this-timestep notch extension cannot exceed the length (i.e. cellside minus notch depth) of sediment remaining on the cell
         pCliff->IncreaseNotchIncision(dNotchIncision);

         // And add to the cell's accumulated wave energy
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->AddToAccumWaveEnergy(dWaveEnergy * dWeight);

         LogStream << m_ulIter << ":\t incision of newly-created notch at [" << nX << "][" << nY << "] dWaveElev = " << dWaveElev << " dCutoffElev = " << dCutoffElev << " dRunup = " << dRunup << "  dWeight = " << dWeight << " dNotchApexElev = " << dNotchApexElev << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dNotchIncision = " << dNotchIncision << endl;

         return true;
      }
      else
      {
         // The top of the notch would be above the top of the sediment on this cell, so we must try to create a notch further inland. Can't do this for points at the beginning and the end of the coast however
         if ((nCoastPoint == 0) || (nCoastPoint == m_VCoast[nCoast].nGetCoastlineSize()-1))
            return false;

         return bCreateNotchInland(nCoast, nCoastPoint, /*nX, nY,*/ dWaveEnergy, dWaveElev);
      }
   }

   // Should never get here
   return false;
}

//===============================================================================================================================
//! If possible, creates an erosional notch further inland from a given coastline point
//===============================================================================================================================
bool CSimulation::bCreateNotchInland(int const nCoast, int const nCoastPoint, /*int const nX, int const nY,*/ double const dWaveEnergy, double const dWaveElev)
{
   // LogStream << m_ulIter << ":\t In bCreateNotchInland() for [" << nX << "][" << nY << "]" << endl;

   bool bFound = false;
   int const nSeaHandedness = m_VCoast[nCoast].nGetSeaHandedness();

   // Get the points before and after this coastline point
   CGeom2DIPoint const* pPtiBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint-1);
   CGeom2DIPoint const* pPtiAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint+1);

   // Loop until we can create an inland notch, or until we hit the edge of the grid
   int n = 1;
   do
   {
      // Get an inland point, planview orthogonal to the coastline at the coast point
      CGeom2DIPoint const PtiTmp = PtiGetPerpendicular(pPtiBefore, pPtiAfter, n, nSeaHandedness);

      // Safety check
      if ((PtiTmp.nGetX() == INT_NODATA) || (! bIsWithinValidGrid(&PtiTmp)))
         return false;

      int const nXTmp = PtiTmp.nGetX();
      int const nYTmp = PtiTmp.nGetY();

      bool bPreExistingNotch;

      // Get the existing notch apex elevation, if there is one
      double dNotchApexElev = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->dGetCliffNotchApexElev();

      if (bFPIsEqual(dNotchApexElev, DBL_NODATA, TOLERANCE))
      {
         // There is no existing notch apex elevation
         bPreExistingNotch = false;

         // So use the notch elevation for this timestep
         dNotchApexElev = m_dThisIterNewNotchApexElev - SED_ELEV_TOLERANCE;
      }
      else
         bPreExistingNotch = true;

      // So, can we create a notch here?
      double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetAllSedTopElevOmitTalus();
      if (dNotchApexElev < dSedTopElevNoTalus)
      {
         // Yes we can potentially create a notch here
         // if (bPreExistingNotch)
         //    LogStream << m_ulIter << ":\t Incision of pre-existing inland cliff notch at [" << nXTmp << "][" << nYTmp << "] dNotchApexElev = " << dNotchApexElev << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << endl;
         // else
         //    LogStream << m_ulIter << ":\t Creation of new notch in inland cliff at [" << nXTmp << "][" << nYTmp << "] dNotchApexElev = " << dNotchApexElev << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << endl;

         // Set the cell to be an inland cliff
         m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->SetLandformCategory(LF_CLIFF);

         // Set its apex elevation
         m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->SetCliffNotchApexElev(dNotchApexElev);

         // We have some notch incision of this newly-created notch
         double dWeight;
         if (dWaveElev > dNotchApexElev)
            // Should not happen: safety check
            dWeight = 1;
         else
            // Assume a linear decrease in incision with distance downwards from notch apex
            dWeight = 1 - ((dNotchApexElev - dWaveElev) / CLIFF_NOTCH_CUTOFF_DISTANCE);

         // Calculate this-timestep cliff notch erosion (is a length in external CRS units). Note that only consolidated sediment can have a cliff notch
         double const dNotchIncision = dWeight * dWaveEnergy / m_dCliffErosionResistance;

         if (bPreExistingNotch)
            m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->AddToCliffNotchIncisionDepth(dNotchIncision);
         else
            m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->SetCliffNotchIncisionDepth(dNotchIncision);

         // And add to the cell's accumulated wave energy
         m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->AddToAccumWaveEnergy(dWaveEnergy * dWeight);

         LogStream << m_ulIter << ":\t inland cliff created (or re-created) at [" << nXTmp << "][" << nYTmp << "] dNotchApexElev = " << dNotchApexElev << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetAllSedTopElevIncTalus() << " increase in notch incision = " << dNotchIncision << " total notch incision = " << m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform()->dGetCliffNotchIncisionDepth() << " threshold incision = " << m_dNotchIncisionAtCollapse << endl;

         bFound = true;
      }

      n++;
   } while (! bFound);

   // Should never get here
   return bFound;
}

//===============================================================================================================================
//! Deposit the unconsolidated sediment from cliff collapse as talus on the cell on which collapse occurred
//===============================================================================================================================
void CSimulation::DoCliffCollapseTalusDeposition(/*int const nCoast,*/ CRWCliff const* pCliff, double const dSandFromCollapse, double const dCoarseFromCollapse, int const nNotchLayer)
{
   // Check: is there some sand- or coarse-sized sediment to deposit?
   if ((dSandFromCollapse + dCoarseFromCollapse) < SED_ELEV_TOLERANCE)
      return;

   // Get the cliff cell's grid coords
   int const nX = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
   int const nY = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();

   // Get a pointer to the layer in which the notch was incised
   CRWCellLayer* pLayer = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer);

   // And get a pointer to the cell layer's talus object
   CRWCellTalus* pTalus = pLayer->pGetOrCreateTalus();

   if (dSandFromCollapse > 0)
   {
      // Add the sand-sized sediment from the collapse to the talus object for this layer
      pTalus->AddSandDepth(dSandFromCollapse);
      m_pRasterGrid->m_Cell[nX][nY].AddSandTalusDeposition(dSandFromCollapse);
   }

   if (dCoarseFromCollapse > 0)
   {
      // Add the coarse-sized sediment from the collapse to the talus object for this layer
      pTalus->AddCoarseDepth(dCoarseFromCollapse);
      m_pRasterGrid->m_Cell[nX][nY].AddCoarseTalusDeposition(dCoarseFromCollapse);
   }

   // And update the cell's sea depth
   m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

   // LogStream << m_ulIter << ":\t coast " << nCoast << " cliff collapse talus deposition on [" << nX << "][" << nY << "] dSandFromCollapse = " << dSandFromCollapse << " dCoarseFromCollapse = " << dCoarseFromCollapse << " sea depth = " << m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() << endl;
}

//===============================================================================================================================
//! Move talus from previous cliff collapse to unconsolidated sediment
//===============================================================================================================================
int CSimulation::nMoveCliffTalusToUnconsolidated(void)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         int const nLayers = m_pRasterGrid->m_Cell[nX][nY].nGetNumLayers();
         for (int nLayer = 0; nLayer < nLayers; nLayer++)
         {
            // Is there talus on this cell layer?
            CRWCellTalus* pTalus = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->pGetTalus();

            if (pTalus == NULL)
               // No talus
               continue;

            // OK we have some talus which could be redistributed from this cell, if waves (inc runup) reach high enough
            double const dThisTalusBottomElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevOmitTalus();

            // Find the nearest point on any coastline, so we can get the runup at this point
            int nCoastClosest;
            int const nCoastPointClosest = nFindClosestCoastPoint(nX, nY, nCoastClosest);
            if (nCoastPointClosest == INT_NODATA)
               return RTN_ERR_CLIFF_TALUS_TO_UNCONS;

            // And get the wave runup at this point
            double const dRunup = m_VCoast[nCoastClosest].dGetRunUp(nCoastPointClosest);

            // Now calc the elevation to which waves reach
            double const dWaveElev = m_dThisIterSWL + dRunup;

            // Only move talus if waves (inc runup) reach above the bottom of the talus
            if (dWaveElev < dThisTalusBottomElev)
            {
               // No talus moved
               if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                  LogStream << m_ulIter << ":\t NO talus moved from [" << nX << "][" << nY << "] since waves do not reach talus base: dWaveElev = " << dWaveElev << " dThisTalusBottomElev = " << dThisTalusBottomElev << endl;

               continue;
            }

            // OK we will move some talus
            double const dThisTalusTopElev = m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus();
            double const dTalusDepth = dThisTalusTopElev - dThisTalusBottomElev;
            double dWeight;

            if (dWaveElev > dThisTalusTopElev)
               // Should not happen: safety check
               dWeight = 1;
            else
               // Assume a linear decrease in talus removal with distance downwards from talus top
               dWeight = 1 - ((dThisTalusTopElev - dWaveElev) / dTalusDepth);

            if (bFPIsEqual(dWeight, 0.0, TOLERANCE))
            {
               // LogStream << m_ulIter << ":\t NO talus moved from [" << nX << "][" << nY << "] dWeight = " << dWeight << endl;

               continue;
            }

            LogStream << m_ulIter << ":\t talus potentially moved from [" << nX << "][" << nY << "] dThisTalusBottomElev = " << dThisTalusBottomElev << " dThisTalusTopElev = " << dThisTalusTopElev << " dWeight = " << dWeight << endl;

            // Next, determine the cells to which talus will be moved. Find all surrounding cells with a top elevation (including talus) which is less than the top elevation (including talus) of this cell
            double dAdjElev;
            double dTotElevDiff = 0;
            vector<double> VdAdjElevDiff;
            vector<CGeom2DIPoint> VptAdj;

            for (int nSearchDirection = NORTH; nSearchDirection <= NORTH_WEST; nSearchDirection++)
            {
               int nXAdj;
               int nYAdj;

               switch (nSearchDirection)
               {
                  case NORTH:
                     nXAdj = nX;
                     nYAdj = nY - 1;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case NORTH_EAST:
                     nXAdj = nX + 1;
                     nYAdj = nY - 1;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case EAST:
                     nXAdj = nX + 1;
                     nYAdj = nY;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case SOUTH_EAST:
                     nXAdj = nX + 1;
                     nYAdj = nY + 1;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case SOUTH:
                     nXAdj = nX;
                     nYAdj = nY + 1;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case SOUTH_WEST:
                     nXAdj = nX - 1;
                     nYAdj = nY + 1;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case WEST:
                     nXAdj = nX - 1;
                     nYAdj = nY;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;

                  case NORTH_WEST:
                     nXAdj = nX - 1;
                     nYAdj = nY - 1;

                     if (bIsWithinValidGrid(nXAdj, nYAdj))
                     {
                        dAdjElev = m_pRasterGrid->m_Cell[nXAdj][nYAdj].dGetAllSedTopElevIncTalus();
                        if (dAdjElev < dThisTalusTopElev)
                        {
                           VptAdj.push_back(CGeom2DIPoint(nXAdj, nYAdj));
                           VdAdjElevDiff.push_back(dAdjElev);
                           dTotElevDiff += dAdjElev;
                        }
                     }

                     break;
               }
            }

            int const nLower = static_cast<int>(VptAdj.size());
            if (nLower == 0)
            {
               // None of the adjacent cells are lower
               // LogStream << m_ulIter << ":\t NO talus moved from [" << nX << "][" << nY << "] since no adjacent cells are lower" << endl;

               continue;
            }

            // OK, at least one adjacent cell is lower. Move talus to each adjacent cell in proportion to the elevation difference
            vector<double> VdPropToMove(nLower);
            for (int n = 0; n < nLower; n++)
               VdPropToMove[n] = VdAdjElevDiff[n] / dTotElevDiff;

            // TODO Removal rate:
            // * to be different for sand and coarse
            // * to include talus erodibility
            // * must depend on SWL and runup.
            // Note we are ignoring subaerial processes.

            double const dTalusSandOrig = pTalus->dGetSandDepth();
            double const dTalusCoarseOrig = pTalus->dGetCoarseDepth();
            double dTalusSandToMove = pTalus->dGetSandDepth();
            double dTalusCoarseToMove = pTalus->dGetCoarseDepth();
            double dTalusSandMoved = 0;
            double dTalusCoarseMoved = 0;
            double const dTalusErodibility = 0.3;
            double const dSandRemovalRate = 1 * dTalusErodibility;            // TEST external CRS units per hour e.g. metres depth per hour (since timestep is in hours)
            double const dCoarseRemovalRate = 0.8 * dTalusErodibility;        // TEST external CRS units per hour e.g. metres depth per hour (since timestep is in hours)

            for (int n = 0; n < nLower; n++)
            {
               if (dTalusSandToMove > 0)
               {
                  // We will deposit some talus sand onto the top layer of this adjacent cell
                  int const nXAdj = VptAdj[n].nGetX();
                  int const nYAdj = VptAdj[n].nGetY();

                  int const nTopLayer = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetNumOfTopLayerAboveBasement();
                  double const dSandOnAdj = m_pRasterGrid->m_Cell[nXAdj][nYAdj].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();

                  double const dPotentialDepthToMove = pTalus->dGetSandDepth() * dWeight * dSandRemovalRate * VdPropToMove[n] * m_dTimeStep;
                  double const dActualDepthToMove = tMin(dTalusSandToMove, dPotentialDepthToMove);

                  m_pRasterGrid->m_Cell[nXAdj][nYAdj].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandOnAdj + dActualDepthToMove);
                  dTalusSandToMove -= dActualDepthToMove;
                  dTalusSandMoved += dActualDepthToMove;

                  // assert(dTalusSandToMove >= 0.0);

                  // Set the changed-this-timestep switch re. the adjacent cell
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  // LogStream << m_ulIter << ":\t " << std::scientific << dActualDepthToMove << std::fixed << " talus sand deposited at [" << nXAdj << "][" << nYAdj << "], talus sand still to deposit on [" << nX << "][" << nY << "] = " << std::scientific << dTalusSandToMove << " talus sand removed = " << dTalusSandMoved << std::fixed << endl;

                  // TODO Update the adjacent cell's talus deposition, and total talus deposition, values
                  // m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dActualDepthToMove);
               }

               if (dTalusCoarseToMove > 0)
               {
                  // We will deposit some talus coarse onto the top layer of this adjacent cell
                  int const nXAdj = VptAdj[n].nGetX();
                  int const nYAdj = VptAdj[n].nGetY();

                  int const nTopLayer = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetNumOfTopLayerAboveBasement();
                  double const dCoarseOnAdj = m_pRasterGrid->m_Cell[nXAdj][nYAdj].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();

                  double const dPotentialDepthToMove = pTalus->dGetCoarseDepth() * dWeight * dCoarseRemovalRate * VdPropToMove[n] * m_dTimeStep;
                  double const dActualDepthToMove = tMin(dTalusCoarseToMove, dPotentialDepthToMove);

                  m_pRasterGrid->m_Cell[nXAdj][nYAdj].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dCoarseOnAdj + dActualDepthToMove);
                  dTalusCoarseToMove -= dActualDepthToMove;
                  dTalusCoarseMoved += dActualDepthToMove;

                  // assert(dTalusCoarseToMove >= 0.0);

                  // LogStream << m_ulIter << ":\t " << std::scientific << dActualDepthToMove << std::fixed << " talus coarse deposited at [" << nXAdj << "][" << nYAdj << "], talus coarse still to deposit on [" << nX << "][" << nY << "] = " << std::scientific << dTalusCoarseToMove << " talus coarse removed = " << dTalusCoarseMoved << std::fixed << endl;

                  // Set the changed-this-timestep switch re. the adjacent cell
                  m_bUnconsChangedThisIter[nTopLayer] = true;

                  // TODO Update the adjacent cell's talus deposition, and total talus deposition, values
                  // m_pRasterGrid->m_Cell[nX][nY].IncrBeachDeposition(dActualDepthToMove);
               }
            }

            if (dTalusSandMoved > 0)
            {
               // For the source cell, update the sand talus value
               double const dTalusSandRemaining = tMax(dTalusSandOrig - dTalusSandMoved, 0.0);

               pTalus->SetSandDepth(dTalusSandRemaining);
            }

            if (dTalusCoarseMoved > 0)
            {
               // For this cell, update the cell layer's coarse talus value
               double const dTalusCoarseRemaining = tMax(dTalusCoarseOrig - dTalusCoarseMoved, 0.0);

               pTalus->SetCoarseDepth(dTalusCoarseRemaining);
            }

            // Has all the talus gone from this layer? If so, then delete it
            if (bFPIsEqual(pTalus->dGetSandDepth() + pTalus->dGetCoarseDepth(), 0.0, TOLERANCE))
            {
               LogStream << m_ulIter << ":\t pTalus->dGetSandDepth() + pTalus->dGetCoarseDepth() = " << std::scientific << pTalus->dGetSandDepth() + pTalus->dGetCoarseDepth() << std::fixed << " so deleting talus object at [" << nX << "][" << nY << "]" << endl;
               m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nLayer)->DeleteTalus();
            }

            // And update the cell's sea depth
            m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

            // LogStream << m_ulIter << ":\t talus moved from [" << nX << "][" << nY << "] sea depth = " << m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() << endl;
         }
      }
   }

   return RTN_OK;
}
