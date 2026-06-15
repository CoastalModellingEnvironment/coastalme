/*!
   \file cell.cpp
   \brief CGeomCell routines
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
   \copyright GNU General Public License
*/

/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <cstddef>

#include <vector>
using std::vector;

#include <assert.h>

#include "cell.h"
#include "cme.h"
#include "simulation.h"
#include "cell.h"
#include "cell_layer.h"
#include "cell_sediment.h"
#include "cell_talus.h"

//! Constructor with initialisation list
CGeomCell::CGeomCell()
   : m_bInContiguousSea(false),
     m_bIsInActiveZone(false),
     // m_bWaveFlood(false),
     // m_bCheckCell(false),
     // m_bCheckFloodCell(false),
     m_bShadowBoundary(false),
     m_bPossibleCoastStartCell(false),
     // m_bPossibleFloodStartCell(false),
     m_bPlatformErosionThisTimestep(false),
     m_nBoundingBoxEdge(NO_DIRECTION),
     m_nCoastlineID(INT_NODATA),
     m_nProfileID(INT_NODATA),
     m_nProfileCoastID(INT_NODATA),
     m_nPolygonID(INT_NODATA),
     m_nPolygonCoastID(INT_NODATA),
     m_nShadowZoneNumber(0),
     m_nDownDriftZoneNumber(0),
     m_dLocalConsSlope(0),
     m_dBasementElevation(0),
     m_dSeaDepth(0),
     m_dTotSeaDepth(0),
     m_dWaveHeight(0),
     m_dTotWaveHeight(0),
     m_dWaveAngle(DBL_NODATA),
     m_dTotWaveAngle(DBL_NODATA),
     m_dDeepWaterWaveHeight(DBL_NODATA),
     m_dDeepWaterWaveAngle(DBL_NODATA),
     m_dBeachProtectionFactor(DBL_NODATA),
     m_dSuspendedSediment(0),
     m_dTotSuspendedSediment(0),
     m_dPotentialPlatformErosionThisIter(0),
     m_dTotPotentialPlatformErosion(0),
     m_dActualPlatformErosionThisIter(0),
     m_dTotActualPlatformErosion(0),
     m_dCliffCollapseToFineTalusThisIter(0),
     m_dCliffCollapseToSandTalusThisIter(0),
     m_dCliffCollapseToCoarseTalusThisIter(0),
     m_dTotFineCliffCollapse(0),
     m_dTotSandCliffCollapse(0),
     m_dTotCoarseCliffCollapse(0),
     m_dTalusFineDepositionThisIter(0),
     m_dTalusSandDepositionThisIter(0),
     m_dTotTalusFineDeposition(0),
     m_dTotTalusSandDeposition(0),
     m_dTalusCoarseDepositionThisIter(0),
     m_dTotTalusCoarseDeposition(0),
     m_dFineTalusToSuspensionThisIter(0),
     m_dTotFineTalusToSuspension(0),
     m_dSandTalusToUnconsThisIter(0),
     m_dTotSandTalusToUncons(0),
     m_dCoarseTalusToUnconsThisIter(0),
     m_dTotCoarseTalusToUncons(0),
     m_dPotentialBeachErosionThisIter(0),
     m_dTotPotentialBeachErosion(0),
     m_dActualBeachErosionThisIter(0),
     m_dTotActualBeachErosion(0),
     m_dBeachDepositionThisIter(0),
     m_dTotBeachDeposition(0),
     m_dUnconsD50(0),
     m_dSlumpingDepositionThisIter(0),
     m_dTotSlumpingDeposition(0),
     m_dInterventionHeight(0)
{
   m_Landform.SetLandformCategory(LF_UNKNOWN);
}

//! Destructor
CGeomCell::~CGeomCell(void)
{
}

//! Set the edge number if this cell is an edge bounding-box cell
void CGeomCell::SetBoundingBoxEdge(int const nDirection)
{
   m_nBoundingBoxEdge = nDirection;
}

//! Returns the number of the bounding-box edge, or NO_DIRECTION if it is not
int CGeomCell::nGetBoundingBoxEdge(void) const
{
   return m_nBoundingBoxEdge;
}

//! Is this an edge bounding-box cell?
bool CGeomCell::bIsBoundingBoxEdge(void) const
{
   return (m_nBoundingBoxEdge != NO_DIRECTION);
}

//! Set this cell to be a sea cell
void CGeomCell::SetInContiguousSea(void)
{
   m_bInContiguousSea = true;
}

//! Is this a sea cell?
bool CGeomCell::bIsInContiguousSea(void) const
{
   return m_bInContiguousSea;
}

// //! Set this cell as flooded by setup surge
// void CGeomCell::SetFloodBySetupSurge(void)
// {
//    m_bFloodBySetupSurge = true;
// }

// //! Is this cell flooded by setup surge?
// bool CGeomCell::bIsFloodBySetupSurge(void) const
// {
//    return m_bFloodBySetupSurge;
// }

// //! Set this cell as flooded by setup surge runup
// void CGeomCell::SetFloodBySetupSurgeRunUp(void)
// {
//    m_bFloodBySetupSurgeRunUp = true;
// }

// //! Is this cell flooded by setup surge runup?
// bool CGeomCell::bIsFloodBySetupSurgeRunUp(void) const
// {
//    return m_bFloodBySetupSurgeRunUp;
// }

//! Sets a flag to show whether this cell is in the active zone
void CGeomCell::SetInActiveZone(bool const bFlag)
{
   m_bIsInActiveZone = bFlag;
}

//! Returns a flag which shows whether this cell is in the active zone
bool CGeomCell::bIsInActiveZone(void) const
{
   return m_bIsInActiveZone;
}

//! Sets a flag to show that this cell is a shadow zone boundary
void CGeomCell::SetShadowZoneBoundary(void)
{
   m_bShadowBoundary = true;
}

//! Returns a flag which shows whether this cell is a shadow zone boundary
bool CGeomCell::bIsShadowZoneBoundary(void) const
{
   return m_bShadowBoundary;
}

//! Sets a flag to show that this cell has been flagged as a possible start- or end-point for a coastline
void CGeomCell::SetPossibleCoastStartCell(void)
{
   m_bPossibleCoastStartCell = true;
}

//! Returns a flag which shows whether this cell has been flagged as a possible start- or end-point for a coastline
bool CGeomCell::bIsPossibleCoastStartCell(void) const
{
   return m_bPossibleCoastStartCell;
}

// //! Sets a flag to show that this cell has been flagged as a possible start-point for runup flooding
// void CGeomCell::SetPossibleFloodStartCell(void)
// {
//    m_bPossibleFloodStartCell = true;
// }

// //! Returns a flag which shows whether this cell has been flagged as a possible start point for runup flooding
// bool CGeomCell::bIsPossibleFloodStartCell(void) const
// {
//    return m_bPossibleFloodStartCell;
// }

//! Returns true if this cell has had potential platform erosion on this cell this timestep
bool CGeomCell::bPotentialPlatformErosion(void) const
{
   return (m_dPotentialPlatformErosionThisIter > 0);
}

//! Returns true if we have had actual platform erosion on this cell during this timestep
bool CGeomCell::bActualPlatformErosion(void) const
{
   return (m_dActualPlatformErosionThisIter > 0);
}

//! Marks this cell with a coastline ID nunber
void CGeomCell::SetAsCoastline(int const nCoast)
{
   m_nCoastlineID = nCoast;
}

//! Returns true if the cell is a coastline
bool CGeomCell::bIsCoastline(void) const
{
   if (m_nCoastlineID == INT_NODATA)
      return false;

   return true;
}

//! Returns the ID of the coastline, or INT_NODATA if not a coastline
int CGeomCell::nGetCoastline(void) const
{
   return m_nCoastlineID;
}

//! Gets the ID number of the coast-normal profile which this cell is 'under', or returns INT_NODATA
int CGeomCell::nGetProfileID(void) const
{
   return m_nProfileID;
}

//! Returns true if this cell is 'under' a coastline normal
bool CGeomCell::bIsProfile(void) const
{
   if (m_nProfileID == INT_NODATA)
      return false;

   return true;
}

//! Gets the coast ID number of the coast-normal profile which this cell is 'under', or returns INT_NODATA
int CGeomCell::nGetProfileCoastID(void) const
{
   return m_nProfileCoastID;
}

//! Sets the coast ID number, and the profile ID number, of the coast-normal profile which this cell is 'under'
void CGeomCell::SetCoastAndProfileID(int const nProfileCoastID, int const nProfile)
{
   m_nProfileCoastID = nProfileCoastID;
   m_nProfileID = nProfile;
}

//! Returns the coast ID number of the polygon which 'contains' this cell (returns INT_NODATA if the cell is not 'in' a polygon)
int CGeomCell::nGetPolygonID(void) const
{
   return m_nPolygonID;
}

//! Returns the coast number of the polygon which 'contains' this cell (returns INT_NODATA if the cell is not 'in' a polygon)
int CGeomCell::nGetPolygonCoastID(void) const
{
   return m_nPolygonCoastID;
}

//! Sets the coast number, and the polygon ID, of the polygon which 'contains' this cell
void CGeomCell::SetCoastAndPolygonID(int const nPolyCoastID, int const nPolyID)
{
   m_nPolygonCoastID = nPolyCoastID;
   m_nPolygonID = nPolyID;
}

//! Set the number of the shadow zone that this cell is in
void CGeomCell::SetShadowZoneNumber(int const nCode)
{
   m_nShadowZoneNumber = nCode;
}

//! Gets the number of the shadow zone that this cell is in
int CGeomCell::nGetShadowZoneNumber(void) const
{
   return m_nShadowZoneNumber;
}

//! Returns true if this cell is in the shadow zone with number given by the parameter, false otherwise
bool CGeomCell::bIsinThisShadowZone(int const nZone) const
{
   if (m_nShadowZoneNumber == nZone)
      return true;

   return false;
}

//! Returns true if this cell is in any shadow zone, false otherwise
bool CGeomCell::bIsinAnyShadowZone(void) const
{
   if (m_nShadowZoneNumber != 0)
      return true;

   return false;
}

//! Returns the depth of talus (all size clases) that is on this cell, in any layer
double CGeomCell::dGetAllTalusDepth(void) const
{
   double dTotTalusDepth = 0;
   for (int nLayer = 0; nLayer < static_cast<int>(m_VLayerAboveBasement.size()); nLayer++)
   {
      CRWCellTalus const* pTalus = m_VLayerAboveBasement[nLayer].pGetTalus();
      if (pTalus != NULL)
      {
         // There is some talus on this layer
         double const dThisTalusDepth = pTalus->dGetFineDepth() + pTalus->dGetSandDepth() + pTalus->dGetCoarseDepth();
         dTotTalusDepth += dThisTalusDepth;
      }
   }

   return dTotTalusDepth;
}

//! Returns the depth of fine talus that is on this cell, in any layer
double CGeomCell::dGetFineTalusDepth(void) const
{
   double dTotTalusDepth = 0;
   for (int nLayer = 0; nLayer < static_cast<int>(m_VLayerAboveBasement.size()); nLayer++)
   {
      CRWCellTalus const* pTalus = m_VLayerAboveBasement[nLayer].pGetTalus();
      if (pTalus != NULL)
      {
         // There is some talus on this layer
         dTotTalusDepth += pTalus->dGetFineDepth();
      }
   }

   return dTotTalusDepth;
}

//! Returns the depth of sand talus that is on this cell, in any layer
double CGeomCell::dGetSandTalusDepth(void) const
{
   double dTotTalusDepth = 0;
   for (int nLayer = 0; nLayer < static_cast<int>(m_VLayerAboveBasement.size()); nLayer++)
   {
      CRWCellTalus const* pTalus = m_VLayerAboveBasement[nLayer].pGetTalus();
      if (pTalus != NULL)
      {
         // There is some talus on this layer
         dTotTalusDepth += pTalus->dGetSandDepth();
      }
   }

   return dTotTalusDepth;
}

//! Returns the depth of coarse talus that is on this cell, in any layer
double CGeomCell::dGetCoarseTalusDepth(void) const
{
   double dTotTalusDepth = 0;
   for (int nLayer = 0; nLayer < static_cast<int>(m_VLayerAboveBasement.size()); nLayer++)
   {
      CRWCellTalus const* pTalus = m_VLayerAboveBasement[nLayer].pGetTalus();
      if (pTalus != NULL)
      {
         // There is some talus on this layer
         dTotTalusDepth += pTalus->dGetCoarseDepth();
      }
   }

   return dTotTalusDepth;
}

// //! Set this cell as flooded by swl + surge + setup + runup
// void CGeomCell::SetWaveFlood(void)
// {
//    m_bWaveFlood = true;
// }

// void CGeomCell::SetWaveSetup(int const dWaveSetup)
// {
// m_dWaveSetup = dWaveSetup;
// }

// void CGeomCell::SetStormSurge(int const dStormSurge)
// {
// m_dStormSurge = dStormSurge;
// }

// void CGeomCell::SetRunUpOnCell(int const dRunUp)
// {
// m_dRunUp = dRunUp;
// }

// void CGeomCell::SetTotLevel(void) const
// {
// m_dTotLevel = m_dWaveSetup + m_dStormSurge + m_dRunUp +
// m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL();
// }

// int CGeomCell::nGetTotLevel(void) const
// {
// return m_dTotLevel;
// }

//! Returns true if the top elevation (sediment and talus, plus any intervention) of this cell is less than this iteration's SWL
bool CGeomCell::bElevLessThanSWL(void) const
{
   // // Note that m_pGrid->pGetSim()->dGetThisIterTotWaterLevel() is zero, since TODO 007 Finish surge and runup stuff
   // return ((m_VdAllHorizonTopElev.back() + this->dGetAllTalusDepth() + m_dInterventionHeight) < (m_pGrid->pGetSim()->dGetThisIterTotWaterLevel() + m_pGrid->pGetSim()->dGetThisIterSWL()));

   // Will need to change this (see above) once surge and runup stuff is working
   return ((m_VdAllHorizonTopElev.back() + this->dGetAllTalusDepth() + m_dInterventionHeight) < m_pGrid->pGetSim()->dGetThisIterSWL());
}

// //! Set this cell as checked TODO 007 Finish surge and runup stuff
// void CGeomCell::SetCheckCell(void)
// {
// m_bCheckCell = true;
// }

// //! Returns true if this cell is checked, false otherwise TODO 007 Finish surge and runup stuff
// bool CGeomCell::bIsCellCheck(void) const
// {
// return m_bCheckCell;
// }

//! Set this cell as checked (flood switch)
void CGeomCell::SetCheckFloodCell(void)
{
   m_bCheckFloodCell = true;
}

//! Set the cell as not checked (flood switch)
void CGeomCell::UnSetCheckFloodCell(void)
{
   m_bCheckFloodCell = false;
}

//! Returns true if this cell is checked, false otherwise (flood switch)
bool CGeomCell::bIsCellFloodCheck(void) const
{
   return m_bCheckFloodCell;
}

//! Sets the down drift zone number
void CGeomCell::SetDownDriftZoneNumber(int const nCode)
{
   m_nDownDriftZoneNumber = nCode;
}

//! Gets the down drift zone number
int CGeomCell::nGetDownDriftZoneNumber(void) const
{
   return m_nDownDriftZoneNumber;
}

//! Returns a pointer to this cell's CRWCellLandform object
CRWCellLandform* CGeomCell::pGetCellLandform(void)
{
   return &m_Landform;
}

//! Sets the local slope of the consolidated sediment only
void CGeomCell::SetLocalConsSedSlope(double const dNewSlope)
{
   m_dLocalConsSlope = dNewSlope;
}

//! Returns the local slope of the consolidated sediment only
double CGeomCell::dGetLocalConsSedSlope(void) const
{
   return m_dLocalConsSlope;
}

//! Sets this cell's basement elevation
void CGeomCell::SetBasementElev(double const dNewElev)
{
   m_dBasementElevation = dNewElev;
}

//! Returns this cell's basement elevation
double CGeomCell::dGetBasementElev(void) const
{
   return (m_dBasementElevation);
}

// //! Sets this cell's slope for cliff toe locating
// void CGeomCell::SetSlopeForCliffToe(double const dNewSlope)
// {
//    m_dSlopeForCliffToe = dNewSlope;
// }
//
// //! Returns this cell's slope for cliff toe locationg
// double CGeomCell::dGetSlopeForCliffToe(void) const
// {
//    return (m_dSlopeForCliffToe);
// }

//! Returns true if this cells's basement data is NODATA, is needed for irregularly-shaped DEMs
bool CGeomCell::bBasementElevIsMissingValue(void) const
{
   if (bFPIsEqual(m_dBasementElevation, m_pGrid->pGetSim()->CSimulation::dGetMissingValue(), TOLERANCE))
      return true;

   return false;
}

//! Returns the depth of seawater on this cell
double CGeomCell::dGetSeaDepth(void) const
{
   return (m_dSeaDepth);
}

//! Returns the total depth of seawater on this cell
double CGeomCell::dGetTotSeaDepth(void) const
{
   return (m_dTotSeaDepth);
}

//! Sets this cell's suspended sediment depth equivalent, it also increments the running total of suspended sediment depth equivalent
void CGeomCell::SetSuspendedSediment(double const dNewSedDepth)
{
   // Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
   m_dSuspendedSediment = dNewSedDepth;
   m_dTotSuspendedSediment += dNewSedDepth;
}

//! Adds to this cell's suspended sediment depth equivalent, it also increments the running total of suspended sediment depth equivalent
void CGeomCell::AddSuspendedSediment(double const dIncSedDepth)
{
   // Note no checks here to see if increment equiv depth is sensible (e.g. non-negative)
   m_dSuspendedSediment += dIncSedDepth;
   m_dTotSuspendedSediment += dIncSedDepth;
}

//! Returns the suspended sediment depth equivalent on this cell
double CGeomCell::dGetSuspendedSediment(void) const
{
   return (m_dSuspendedSediment);
}

//! Returns the total suspended sediment depth equivalent on this cell
double CGeomCell::dGetTotSuspendedSediment(void) const
{
   return (m_dTotSuspendedSediment);
}

//! Returns the number of sediment layers on the cell (some could have zero thickness)
int CGeomCell::nGetNumLayers(void) const
{
   return static_cast<int>(m_VLayerAboveBasement.size());
}


//! Returns the index of the topmost sediment layer (layer 0 being the one just above basement) with non-zero thickness. If there is no such layer, it returns NO_NONZERO_THICKNESS_LAYERS
int CGeomCell::nGetTopNonZeroLayerAboveBasement(void) const
{
   if (m_VLayerAboveBasement.empty())
      return INT_NODATA;

   int nTop = static_cast<int>(m_VLayerAboveBasement.size()) - 1;

   while (m_VLayerAboveBasement[nTop].dGetTotalThickness() <= 0)
   {
      if (--nTop < 0)
         return NO_NONZERO_THICKNESS_LAYERS;
   }

   return nTop;
}

//! Returns the index of the topmost sediment layer (layer 0 being the one just above basement). Note that this layer could have zero thickness
int CGeomCell::nGetNumOfTopLayerAboveBasement(void) const
{
   if (m_VLayerAboveBasement.empty())
      return INT_NODATA;

   return static_cast<int>(m_VLayerAboveBasement.size()) - 1;
}

//! Returns the elevation of the top of the consolidated sediment only, for a given layer (layer 0 being the one just above basement)
double CGeomCell::dGetConsSedTopElevForLayerAboveBasement(int const nLayer) const
{
   // Note no check to see if nLayer < m_VLayerAboveBasement.size()
   double dTopElev = m_dBasementElevation;

   for (int n = 0; n < nLayer; n++)
   {
      dTopElev += m_VLayerAboveBasement[n].dGetAllUnconsDepth();
      dTopElev += m_VLayerAboveBasement[n].dGetAllConsDepth();
   }

   dTopElev += m_VLayerAboveBasement[nLayer].dGetAllConsDepth();

   return dTopElev;
}

//! Return a reference to the Nth sediment layer (layer 0 being just above basement)
CRWCellLayer* CGeomCell::pGetLayerAboveBasement(int const nLayer)
{
   // TODO 055 No check that nLayer < size()
   return &m_VLayerAboveBasement[nLayer];
}

//! Returns the elevation of the top surface of sediment (both consolidated and unconsolidated) for this cell, but ignoring any talus. If there is a cliff notch, ignore the missing volume
double CGeomCell::dGetAllSedTopElevOmitTalus(void) const
{
   return m_VdAllHorizonTopElev.back();
}

//! Returns the elevation of the top surface of sediment (both consolidated and unconsolidated) plus talus for this cell. If there is a cliff notch, ignore the missing volume
double CGeomCell::dGetAllSedTopElevIncTalus(void)
{
   return this->dGetAllSedTopElevOmitTalus() + this->dGetAllTalusDepth();
}

//! Returns the elevation of the top surface of sediment (both consolidated and unconsolidated) for this cell, but ignoring any talus. If there is a cliff notch, ignore the missing volume
double CGeomCell::dGetConsSedTopElevOmitTalus(void) const
{
   return this->dGetBasementElev() + this->dGetAllConsDepthAllLayers();
}

//! Returns the elevation of the top surface of sediment (both consolidated and unconsolidated) plus talus for this cell. If there is a cliff notch, ignore the missing volume
double CGeomCell::dGetConsSedTopElevIncTalus(void)
{
   return this->dGetConsSedTopElevOmitTalus() + this->dGetAllTalusDepth();
}

//! Returns the topmost elevation of the cell, including sea (sediment top elevation (both consolidated and unconsolidated), plus the depth of any talus, plus the height of any intervention, plus the sea depth)
double CGeomCell::dGetTopElevIncSea(void)
{
   return m_VdAllHorizonTopElev.back() + this->dGetAllTalusDepth() + m_dInterventionHeight + m_dSeaDepth;
}

//! Returns true if the elevation of the sediment top surface (both consolidated and unconsolidated, and any talus) for this cell, plus any intervention, is less than the grid's this-timestep still water elevation
bool CGeomCell::bIsInundated(void)
{
   return ((m_VdAllHorizonTopElev.back() + m_dInterventionHeight + this->dGetAllTalusDepth()) < m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL());
}

//! Returns the sea surface elevation at the current iteration
double CGeomCell::dGetThisIterSWL(void) const
{
   return m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL();
}

//! Returns the total water level for the current iteration  TODO 007 Finish surge and runup stuff
double CGeomCell::dGetThisIterTotWaterLevel(void) const
{
   return m_pGrid->pGetSim()->CSimulation::dGetThisIterTotWaterLevel();
}

//! Returns the total (all layers) depth (in external CRS units) of fine consolidated sediment on this cell
double CGeomCell::dGetConsFineDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetConsFineDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of fine unconsolidated sediment on this cell
double CGeomCell::dGetUnconsFineDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetUnconsFineDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of sand-sized consolidated sediment on this cell
double CGeomCell::dGetConsSandDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetConsSandDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of sand-sized unconsolidated sediment on this cell
double CGeomCell::dGetUnconsSandDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetUnconsSandDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of coarse consolidated sediment on this cell
double CGeomCell::dGetConsCoarseDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetConsCoarseDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of coarse unconsolidated sediment on this cell
double CGeomCell::dGetUnconsCoarseDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetUnconsCoarseDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of consolidated sediment (all size classes) on this cell
double CGeomCell::dGetAllConsDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetAllConsDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of unconsolidated sediment (all size classes) on this cell
double CGeomCell::dGetAllUnconsDepthAllLayers(void) const
{
   double dTotDepth = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotDepth += m_VLayerAboveBasement[n].dGetAllUnconsDepth();

   return dTotDepth;
}

//! Returns the total (all layers) depth (in external CRS units) of all sediment (all size classes, both consolidated and unconsolidated) on this cell
double CGeomCell::dGetAllSedDepthAllLayers(void) const
{
   return (this->dGetAllUnconsDepthAllLayers() + this->dGetAllConsDepthAllLayers());
}

//! Appends sediment layers
void CGeomCell::AppendLayers(int const nLayer)
{
   for (int i = 0; i < nLayer; i++)
      m_VLayerAboveBasement.push_back(CRWCellLayer());
}

//! For this cell, calculates the elevation of the top of every layer, and the d50 for the topmost unconsolidated sediment layer
void CGeomCell::CalcAllLayerElevsAndD50(void)
{
   m_VdAllHorizonTopElev.clear();
   m_VdAllHorizonTopElev.push_back(m_dBasementElevation);      // Elevation of top of the basement

   // Calculate the elevation of the top of all other layers
   int m = 0;

   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      m_VdAllHorizonTopElev.push_back(m_VLayerAboveBasement[n].dGetTotalThickness() + m_VdAllHorizonTopElev[m++]); // Elevation of top of layer n, ignoring any talus

   // Now calculate the d50 of the topmost unconsolidated sediment layer with non-zero thickness. If there is no unconsolidated sediment, m_dUnconsD50 is set to DBL_NODATA
   m_dUnconsD50 = DBL_NODATA;

   for (int n = static_cast<int>(m_VLayerAboveBasement.size()) - 1; n >= 0; n--)
   {
      double const dUnconsThick = m_VLayerAboveBasement[n].dGetAllUnconsDepth();

      if (dUnconsThick > 0)
      {
         // This is a layer with non-zero thickness of unconsolidated sediment
         CRWCellSediment const* pUnconsSedLayer = m_VLayerAboveBasement[n].pGetUnconsolidatedSediment();
         double const dFineProp = pUnconsSedLayer->dGetFineDepth() / dUnconsThick;
         double const dSandProp = pUnconsSedLayer->dGetSandDepth() / dUnconsThick;
         double const dCoarseProp = pUnconsSedLayer->dGetCoarseDepth() / dUnconsThick;

         // Calculate d50 for the unconsolidated sediment
         m_dUnconsD50 = (dFineProp * m_pGrid->pGetSim()->dGetD50Fine()) + (dSandProp * m_pGrid->pGetSim()->dGetD50Sand()) + (dCoarseProp * m_pGrid->pGetSim()->dGetD50Coarse());

         break;
      }
   }
}

//! Given an elevation, this finds the index of the layer that contains that elevation (layer 0 being the first above the basement). Note that the elevation cannot exactly equal the elevation of the layer's top surface
int CGeomCell::nGetLayerAtElev(double const dElev) const
{
   /*! Returns ELEV_IN_BASEMENT if in basement, ELEV_ABOVE_SEDIMENT_TOP if higher than or equal to sediment top, or the layer number (0 to n),  */
   if (dElev < m_VdAllHorizonTopElev[0])
      return ELEV_IN_BASEMENT;

   for (unsigned int nLayer = 1; nLayer < m_VdAllHorizonTopElev.size();
        nLayer++)
   {
      if ((m_VLayerAboveBasement[nLayer - 1].dGetTotalThickness() > 0) && (dElev >= m_VdAllHorizonTopElev[nLayer - 1]) && (dElev <= m_VdAllHorizonTopElev[nLayer]))
         return (nLayer - 1);
   }

   return ELEV_ABOVE_SEDIMENT_TOP;
}

//! For this cell, calculates the elevation of the top of a given layer
double CGeomCell::dCalcLayerElev(const int nLayer)
{
   // Note no check to see if nLayer < m_VLayerAboveBasement.size()
   double dTopElev = m_dBasementElevation;

   for (int n = 0; n <= nLayer; n++)
      dTopElev += m_VLayerAboveBasement[n].dGetTotalThickness();

   return dTopElev;
}

//! Set potential (unconstrained) shore platform erosion and increment total shore platform potential erosion
void CGeomCell::SetPotentialPlatformErosion(double const dPotentialIn)
{
   m_dPotentialPlatformErosionThisIter = dPotentialIn;
   m_dTotPotentialPlatformErosion += dPotentialIn;
}

//! Get potential (unconstrained) shore platform erosion
double CGeomCell::dGetPotentialPlatformErosion(void) const
{
   return m_dPotentialPlatformErosionThisIter;
}

//! Get total potential (unconstrained) shore platform erosion
double CGeomCell::dGetTotPotentialPlatformErosion(void) const
{
   return m_dTotPotentialPlatformErosion;
}

//! Set this-timestep actual (constrained) shore platform erosion and increment total actual shore platform erosion
void CGeomCell::SetActualPlatformErosion(double const dThisActualErosion)
{
   m_dActualPlatformErosionThisIter = dThisActualErosion;
   m_dTotActualPlatformErosion += dThisActualErosion;
}

//! Get actual (constrained) shore platform erosion
double CGeomCell::dGetActualPlatformErosion(void) const
{
   return m_dActualPlatformErosionThisIter;
}

//! Get total actual (constrained) shore platform erosion
double CGeomCell::dGetTotActualPlatformErosion(void) const
{
   return m_dTotActualPlatformErosion;
}

//! Returns the depth of seawater on this cell if the sediment top (including talus and intervention) is < SWL, or zero
void CGeomCell::SetSeaDepth(void)
{
   m_dSeaDepth = tMax(m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL() - (m_VdAllHorizonTopElev.back() + this->dGetAllTalusDepth() + m_dInterventionHeight), 0.0);
}

//! Initialise values for this cell
void CGeomCell::InitCell(void)
{
   m_bInContiguousSea = false;
   m_bIsInActiveZone = false;
   // m_bEstimated = false;
   m_bShadowBoundary = false;
   m_bPossibleCoastStartCell = false;
   // m_bPossibleFloodStartCell = false;
   // m_bWaveFlood = false;
   // m_bCheckCell = false;
   // m_bCheckFloodCell = false;

   m_nCoastlineID = INT_NODATA;
   m_nPolygonID = INT_NODATA;
   m_nPolygonCoastID = INT_NODATA;
   m_nProfileID = INT_NODATA;

   m_nShadowZoneNumber = 0;
   m_nDownDriftZoneNumber = 0;

   m_dLocalConsSlope = 0;
   m_dPotentialPlatformErosionThisIter = 0;
   m_dActualPlatformErosionThisIter = 0;
   m_dCliffCollapseToFineTalusThisIter = 0;
   m_dCliffCollapseToSandTalusThisIter = 0;
   m_dCliffCollapseToCoarseTalusThisIter = 0;
   m_dTalusFineDepositionThisIter = 0;
   m_dTalusSandDepositionThisIter = 0;
   m_dTalusCoarseDepositionThisIter = 0;
   m_dFineTalusToSuspensionThisIter = 0;
   m_dSandTalusToUnconsThisIter = 0;
   m_dCoarseTalusToUnconsThisIter = 0;
   m_dPotentialBeachErosionThisIter = 0;
   m_dActualBeachErosionThisIter = 0;
   m_dBeachDepositionThisIter = 0;
   m_dSeaDepth = 0;
   m_dWaveHeight = 0;
   m_dWaveAngle = 0;

   m_dBeachProtectionFactor = DBL_NODATA;

   // Initialize this-iteration sediment input event values
   int const nThisLayer = this->nGetTopNonZeroLayerAboveBasement();

   // Safety check
   if ((nThisLayer == NO_NONZERO_THICKNESS_LAYERS) ||
       (nThisLayer == INT_NODATA))
   {
      // Uh-oh, problem. So just return TODO 021
      return;
   }

   m_VLayerAboveBasement[nThisLayer].pGetUnconsolidatedSediment()->InitThisIterSedimentInputAll();
}

//! Sets the wave height on this cell
void CGeomCell::SetWaveHeight(double const dWaveHeight)
{
   // assert(m_dWaveHeight != DBL_NODATA);
   // assert(m_dWaveHeight >= 0);
   m_dWaveHeight = dWaveHeight;
}

//! Returns the wave height on this cell
double CGeomCell::dGetWaveHeight(void) const
{
   return m_dWaveHeight;
}

//! Returns the total wave height on this cell
double CGeomCell::dGetTotWaveHeight(void) const
{
   return m_dTotWaveHeight;
}

//! Increments the total wave height on this cell
void CGeomCell::IncrTotWaveHeight(const double dWaveHeight)
{
   m_dTotWaveHeight += dWaveHeight;
}

//! Sets the wave orientation on this cell
void CGeomCell::SetWaveAngle(double const dWaveAngle)
{
   m_dWaveAngle = dWaveAngle;
}

//! Returns the wave orientation on this cell
double CGeomCell::dGetWaveAngle(void) const
{
   return m_dWaveAngle;
}

//! Increments to total wave orientation on this cell
void CGeomCell::IncrTotWaveAngle(double const dWaveAngle)
{
   m_dTotWaveAngle += dWaveAngle;
}

//! Returns the total wave orientation on this cell
double CGeomCell::dGetTotWaveAngle(void) const
{
   return m_dTotWaveAngle;
}

//! Sets the deep water wave height on this cell
void CGeomCell::SetCellDeepWaterWaveHeight(double const dWaveHeight)
{
   m_dDeepWaterWaveHeight = dWaveHeight;
}

//! Returns the deep water wave height on this cell
double CGeomCell::dGetCellDeepWaterWaveHeight(void) const
{
   return m_dDeepWaterWaveHeight;
}

//! Sets the deep water wave orientation on this cell
void CGeomCell::SetCellDeepWaterWaveAngle(double const dWaveAngle)
{
   m_dDeepWaterWaveAngle = dWaveAngle;
}

//! Returns the deep water wave orientation on this cell
double CGeomCell::dGetCellDeepWaterWaveAngle(void) const
{
   return m_dDeepWaterWaveAngle;
}

//! Sets the deep water wave Period on this cell
void CGeomCell::SetCellDeepWaterWavePeriod(double const dWavePeriod)
{
   m_dDeepWaterWavePeriod = dWavePeriod;
}

//! Returns the deep water wave period on this cell
double CGeomCell::dGetCellDeepWaterWavePeriod(void) const
{
   return m_dDeepWaterWavePeriod;
}

//! Sets wave height to the deep water wave height value, and sets wave orientation to the deep water wave orientation value
void CGeomCell::SetWaveValuesToDeepWaterWaveValues(void)
{
   m_dWaveHeight = m_dDeepWaterWaveHeight;
   m_dWaveAngle = m_dDeepWaterWaveAngle;
   m_dWavePeriod = m_dDeepWaterWavePeriod;
}

// Sets this cell's beach protection factor
void CGeomCell::SetBeachProtectionFactor(double const dFactor)
{
   m_dBeachProtectionFactor = dFactor;
}

//! Returns this cell's beach protection factor
double CGeomCell::dGetBeachProtectionFactor(void) const
{
   return m_dBeachProtectionFactor;
}

//! Increments the fine, sand, and coarse depths of this-timestep cliff collapse on this cell, also increments the totals
void CGeomCell::IncrCliffCollapseErosion(double const dFineDepth, double const dSandDepth, double const dCoarseDepth)
{
   m_dCliffCollapseToFineTalusThisIter += dFineDepth;
   m_dCliffCollapseToSandTalusThisIter += dSandDepth;
   m_dCliffCollapseToCoarseTalusThisIter += dCoarseDepth;

   m_dTotFineCliffCollapse += dFineDepth;
   m_dTotSandCliffCollapse += dSandDepth;
   m_dTotCoarseCliffCollapse += dCoarseDepth;
}

//! Returns the depth of this-timestep cliff collapse fine-sized sediment moved to talus on this cell
double CGeomCell::dGetThisIterCliffCollapseErosionFineToTalus(void) const
{
   return m_dCliffCollapseToFineTalusThisIter;
}

//! Returns the depth of this-timestep cliff collapse sand-sized sediment moved to talus on this cell
double CGeomCell::dGetThisIterCliffCollapseErosionSandToTalus(void) const
{
   return m_dCliffCollapseToSandTalusThisIter;
}

//! Returns the depth of this-timestep cliff collapse coarse-sized sediment moved to talus on this cell
double CGeomCell::dGetThisIterCliffCollapseErosionCoarseToTalus(void) const
{
   return m_dCliffCollapseToCoarseTalusThisIter;
}

//! Returns the running total depth of fine-sized sediment eroded by cliff collapse on this cell
double CGeomCell::dGetTotCliffCollapseFine(void) const
{
   return m_dTotFineCliffCollapse;
}

//! Returns the running total depth of sand-sized sediment eroded by cliff collapse on this cell
double CGeomCell::dGetTotCliffCollapseSand(void) const
{
   return m_dTotSandCliffCollapse;
}

//! Returns the running total depth of coarse-sized sediment eroded by cliff collapse on this cell
double CGeomCell::dGetTotCliffCollapseCoarse(void) const
{
   return m_dTotCoarseCliffCollapse;
}

//! Increments the depth of this-timestep fine-sized cliff collapse talus on this cell, also increments the total
void CGeomCell::AddFineTalusDeposition(double const dDepth)
{
   m_dTalusFineDepositionThisIter += dDepth;
   m_dTotTalusFineDeposition += dDepth;
}

//! Increments the depth of this-timestep sand-sized cliff collapse talus on this cell, also increments the total
void CGeomCell::AddSandTalusDeposition(double const dDepth)
{
   m_dTalusSandDepositionThisIter += dDepth;
   m_dTotTalusSandDeposition += dDepth;
}

//! Increments the depth of this-timestep coarse-sized cliff collapse talus on this cell, also increments the total
void CGeomCell::AddCoarseTalusDeposition(double const dDepth)
{
   m_dTalusCoarseDepositionThisIter += dDepth;
   m_dTotTalusCoarseDeposition += dDepth;
}

//! Returns the depth of this-timestep sand talus deposition from cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseSandTalusDeposition(void) const
{
   return m_dTalusSandDepositionThisIter;
}

//! Returns the depth of this-timestep coarse talus deposition from cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseCoarseTalusDeposition(void) const
{
   return m_dTalusCoarseDepositionThisIter;
}

//! Returns the total depth of sand talus deposition from cliff collapse on this cell
double CGeomCell::dGetTotSandTalusDeposition(void) const
{
   return m_dTotTalusSandDeposition;
}

//! Returns the total depth of coarse talus deposition from cliff collapse on this cell
double CGeomCell::dGetTotCoarseTalusDeposition(void) const
{
   return m_dTotTalusCoarseDeposition;
}

//! Adds to the totals of fine-sized talus moved to suspension on this cell
void CGeomCell::AddFineTalusToSuspension(double const dFineTalus)
{
   m_dFineTalusToSuspensionThisIter += dFineTalus;
   m_dTotFineTalusToSuspension += dFineTalus;

}

//! Returns the this-iteration total of fine-sized talus moved to suspension on this cell
double CGeomCell::dGetThisIterFineTalusToSuspension(void)
{
   return m_dFineTalusToSuspensionThisIter;
}

//! Returns the total depth of fine-sized talus moved to suspension on this cell
double CGeomCell::dGetTotFineTalusToSuspension(void)
{
   return m_dTotFineTalusToSuspension;
}

//! Adds to the totals of sand-sized talus moved to unconsolidated sediment on this cell
void CGeomCell::AddSandTalusToUncons(double const dSandTalus)
{
   m_dSandTalusToUnconsThisIter += dSandTalus;
   m_dTotSandTalusToUncons += dSandTalus;
}

//! Returns the this-iteration total of sand-sized talus moved to unconsolidated sediment on this cell
double CGeomCell::dGetThisIterSandTalusToUncons(void)
{
   return m_dSandTalusToUnconsThisIter;
}

//! Returns the grand total of sand-sized talus moved to unconsolidated sediment on this cell
double CGeomCell::dGetTotSandTalusToUncons(void)
{
   return m_dTotSandTalusToUncons;
}

//! Adds to the totals of coarse-sized talus moved to unconsolidated sediment on this cell
void CGeomCell::AddCoarseTalusToUncons(double const dCoarseTalus)
{
   m_dCoarseTalusToUnconsThisIter += dCoarseTalus;
   m_dTotCoarseTalusToUncons += dCoarseTalus;
}

//! Returns the this-iteration total of coarse-sized talus moved to unconsolidated sediment on this cell
double CGeomCell::dGetThisIterCoarseTalusToUncons(void)
{
   return m_dCoarseTalusToUnconsThisIter;
}

//! Returns the grand total of coarse-sized talus moved to unconsolidated sediment on this cell
double CGeomCell::dGetTotCoarseTalusToUncons(void)
{
   return m_dTotCoarseTalusToUncons;
}

//! Set potential (unconstrained) beach erosion and increment total beach potential erosion
void CGeomCell::SetPotentialBeachErosion(double const dPotentialIn)
{
   m_dPotentialBeachErosionThisIter = dPotentialIn;
   m_dTotPotentialBeachErosion += dPotentialIn;
}

//! Get potential (unconstrained) beach erosion
double CGeomCell::dGetPotentialBeachErosion(void) const
{
   return m_dPotentialBeachErosionThisIter;
}

//! Get total potential (supply-unconstrained) beach erosion
double CGeomCell::dGetTotPotentialBeachErosion(void) const
{
   return m_dTotPotentialBeachErosion;
}

//! Set this-timestep actual (supply-constrained) beach erosion and increment total actual beach erosion
void CGeomCell::SetActualBeachErosion(double const dThisActualErosion)
{
   m_dActualBeachErosionThisIter = dThisActualErosion;
   m_dTotActualBeachErosion += dThisActualErosion;
}

//! Get actual (supply-constrained) beach erosion
double CGeomCell::dGetActualBeachErosion(void) const
{
   return m_dActualBeachErosionThisIter;
}

//! Get total actual (supply-constrained) beach erosion
double CGeomCell::dGetTotActualBeachErosion(void) const
{
   return m_dTotActualBeachErosion;
}

// //! Returns true if there has been actual beach erosion this timestep
// bool CGeomCell::bActualBeachErosionThisIter(void) const
// {
// return (m_dActualBeachErosionThisIter > 0 ? true : false);
// }

//! Increment this-timestep beach deposition, also increment total beach deposition
void CGeomCell::IncrBeachDeposition(double const dThisDeposition)
{
   m_dBeachDepositionThisIter += dThisDeposition;
   m_dTotBeachDeposition += dThisDeposition;
}

//! Get beach deposition
double CGeomCell::dGetBeachDeposition(void) const
{
   return m_dBeachDepositionThisIter;
}

//! Get beach erosion
double CGeomCell::dGetTotBeachDeposition(void) const
{
   return m_dTotBeachDeposition;
}

// //! Returns true if there has been beach deposition this timestep
// bool CGeomCell::bBeachDepositionThisIter(void) const
// {
// return (m_dBeachDepositionThisIter > 0 ? true : false);
// }

//! Returns true only if this cell has had no deposition or erosion this timestep
bool CGeomCell::bBeachErosionOrDepositionThisIter(void) const
{
   if ((m_dActualBeachErosionThisIter > 0) || (m_dBeachDepositionThisIter > 0))
      return true;

   return false;
}

//! Returns the D50 of unconsolidated sediment on this cell
double CGeomCell::dGetUnconsD50(void) const
{
   return m_dUnconsD50;
}

//! Increment this-timestep avalanche deposition, also increment total avalanche deposition
void CGeomCell::IncrSlumpDeposition(double const dThisDeposition)
{
   m_dSlumpingDepositionThisIter += dThisDeposition;
   m_dTotSlumpingDeposition += dThisDeposition;
}

//! Get avalanche deposition for this timestep
double CGeomCell::dGetAvalancheDeposition(void) const
{
   return m_dSlumpingDepositionThisIter;
}

//! Get total avalanche deposition
double CGeomCell::dGetTotAvalancheDeposition(void) const
{
   return m_dTotSlumpingDeposition;
}

//! Gets the intervention class for this cell, or returns INT_NODATA if there is not an intervention here
int CGeomCell::nGetInterventionClass(void) const
{
   int const nTmp = this->m_Landform.nGetLandformCategory();
   if ((nTmp == LF_INTERVENTION_STRUCT) || (nTmp == LF_INTERVENTION_NON_STRUCT))
      return nTmp;

   return INT_NODATA;
}

//! Sets the intervention height
void CGeomCell::SetInterventionHeight(double const dHeight)
{
   m_dInterventionHeight = dHeight;
}

//! Returns the intervention height
double CGeomCell::dGetInterventionHeight(void) const
{
   return m_dInterventionHeight;
}

//! Returns the elevation of the top of the intervention, assuming it rests on the sediment-top surface
double CGeomCell::dGetInterventionTopElev(void) const
{
   return m_VdAllHorizonTopElev.back() + m_dInterventionHeight;
}

//! Set the switch to show platform erosion this timestep
void CGeomCell::SetPlatformErosionThisIter(bool const bSwitch)
{
   m_bPlatformErosionThisTimestep = bSwitch;
}

//! Returns the value of the switch that shows platform erosion this timestep
bool CGeomCell::bGetPlatformErosionThisIter(void) const
{
   return m_bPlatformErosionThisTimestep;
}
