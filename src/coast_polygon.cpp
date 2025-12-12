/*!
   \file coast_polygon.cpp
   \brief CGeomCoastPolygon routines
   \details TODO 001 A more detailed description of these routines.
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

#include <vector>
using std::vector;

#include <algorithm>
using std::find;

#include "cme.h"
#include "coast_polygon.h"
#include "2d_point.h"
#include "2di_point.h"

//! Constructor with 10 parameters and initialisation list
CGeomCoastPolygon::CGeomCoastPolygon(int const nCoastID, int const nNode, int const nProfileUpCoast, int const nProfileDownCoast, vector<CGeom2DPoint> const* pVIn, int const nLastPointUpCoast, const int nLastPointDownCoast, CGeom2DIPoint const* PtiNode, CGeom2DIPoint const* PtiAntinode, bool const bStartCoast, bool const bEndCoast)
    : // m_bIsPointedSeaward(true),
      m_bUnconsSedimentMovementDownCoastThisIter(false),
      m_bCoastEndPolygon(bEndCoast),
      m_bCoastStartPolygon(bStartCoast),
      m_nCoastID(nCoastID),
      m_nCoastNode(nNode),
      m_nProfileUpCoast(nProfileUpCoast),
      m_nProfileDownCoast(nProfileDownCoast),
      m_nProfileUpCoastNumPointsUsed(nLastPointUpCoast),
      m_nProfileDownCoastNumPointsUsed(nLastPointDownCoast),
      m_nNumCells(0),
      m_dAvgUnconsD50(0),
      m_dSeawaterVolume(0),
      m_dPotentialBeachErosionAllUncons(0),
      m_dBeachErosionUnconsFine(0),
      m_dBeachErosionUnconsSand(0),
      m_dBeachErosionUnconsCoarse(0),
      m_dBeachDepositionUnconsSand(0),
      m_dBeachDepositionUnconsCoarse(0),
      m_dSuspensionUnconsFine(0),
      m_dToDoBeachDepositionUnconsSand(0),
      m_dToDoBeachDepositionUnconsCoarse(0),
      m_dBeachSandErodedDeanProfile(0),
      m_dBeachCoarseErodedDeanProfile(0),
      m_dCliffCollapseErosionFine(0),
      m_dCliffCollapseErosionSand(0),
      m_dCliffCollapseErosionCoarse(0),
      m_dCliffCollapseTalusFineToSuspension(0),
      m_dCliffCollapseTalusSand(0),
      m_dCliffCollapseTalusCoarse(0),
      m_dCliffCollapseSandErodedDeanProfile(0),
      m_dCliffCollapseCoarseErodedDeanProfile(0),
      m_dPlatformErosionToSuspensionFine(0),
      m_dPlatformErosionUnconsSand(0),
      m_dPlatformErosionUnconsCoarse(0),
      m_dPreExistingUnconsFine(0),
      m_dPreExistingUnconsSand(0),
      m_dPreExistingUnconsCoarse(0),
      m_dPreExistingConsFine(0),
      m_dPreExistingConsSand(0),
      m_dPreExistingConsCoarse(0),
      m_dSedimentInputFine(0),
      m_dSedimentInputSand(0),
      m_dSedimentInputCoarse(0),
      m_dLength(0),
      m_PtiNode(*PtiNode),
      m_PtiAntinode(*PtiAntinode)
{
   m_VPoints = *pVIn;
}

//! Destructor
CGeomCoastPolygon::~CGeomCoastPolygon(void)
{
}

// void CGeomCoastPolygon::SetNotPointed(void)
// {
// m_bIsPointedSeaward = false;
// }
//
// bool CGeomCoastPolygon::bIsPointed(void) const
// {
// return m_bIsPointedSeaward;
// }

//! Set a flag to say whether sediment movement on this polygon is down-coast this iteration
void CGeomCoastPolygon::SetDownCoastThisIter(bool const bFlag)
{
   m_bUnconsSedimentMovementDownCoastThisIter = bFlag;
}

//! Is sediment movement on this polygon down-coast this iteration?
bool CGeomCoastPolygon::bDownCoastThisIter(void) const
{
   return m_bUnconsSedimentMovementDownCoastThisIter;
}

// //! Set this coast polygon as the coast-end polygon
// void CGeomCoastPolygon::SetCoastEndPolygon(void)
// {
//    m_bCoastEndPolygon = true;
// }

//! Is this polygon the coast-end polygon?
bool CGeomCoastPolygon::bIsCoastEndPolygon(void) const
{
   return m_bCoastEndPolygon;
}

// //! Set this coast polygon as the coast-start polygon
// void CGeomCoastPolygon::SetCoastStartPolygon(void)
// {
//    m_bCoastStartPolygon = true;
// }

//! Is this polygon the coast-start polygon?
bool CGeomCoastPolygon::bIsCoastStartPolygon(void) const
{
   return m_bCoastStartPolygon;
}

//! Get the coast ID, this is the same as the down-coast sequence of polygons
int CGeomCoastPolygon::nGetPolygonCoastID(void) const
{
   return m_nCoastID;
}

// void CGeomCoastPolygon::SetCoastNode(int const nNode)
// {
// m_nCoastNode = nNode;
// }

//! Get the coast node point
int CGeomCoastPolygon::nGetNodeCoastPoint(void) const
{
   return m_nCoastNode;
}

//! Get the grid coordinates of the cell on which the node sits
CGeom2DIPoint* CGeomCoastPolygon::pPtiGetNode(void)
{
   return &m_PtiNode;
}

//! Get the anti-node (raster grid CRS) which is at other (seaward) end of the polygon from the node
CGeom2DIPoint* CGeomCoastPolygon::pPtiGetAntiNode(void)
{
   return &m_PtiAntinode;
}

//! Sets the polygon's length
void CGeomCoastPolygon::SetLength(double const dLength)
{
   m_dLength = dLength;
}

//! Gets the polygon's length
double CGeomCoastPolygon::dGetLength(void) const
{
   return m_dLength;
}

//! Set the number of cells in the polygon
void CGeomCoastPolygon::SetNumCellsInPolygon(int const nCells)
{
   m_nNumCells = nCells;
}

// //! Get the number of cells in the polygon
// int CGeomCoastPolygon::nGetNumCellsinPolygon(void) const
// {
// return m_nNumCells;
// }

//! Return the number of the up-coast profile
int CGeomCoastPolygon::nGetUpCoastProfile(void) const
{
   return m_nProfileUpCoast;
}

//! Return the number of the down-coast profile
int CGeomCoastPolygon::nGetDownCoastProfile(void) const
{
   return m_nProfileDownCoast;
}

// void CGeomCoastPolygon::SetBoundary(vector<CGeom2DPoint> const* pVIn)
// {
// m_VPoints = *pVIn;
// }

// vector<CGeom2DPoint>* CGeomCoastPolygon::pPtVGetBoundary(void)
// {
// return &m_VPoints;
// }

//! Get the coordinates (external CRS) of a specified point on the polygon's boundary
CGeom2DPoint* CGeomCoastPolygon::pPtGetBoundaryPoint(int const nPoint)
{
   // TODO 055 No check to see if nPoint < m_VPoints.size()
   return &m_VPoints[nPoint];
}

//! Get the number of points in the polygon's boundary
int CGeomCoastPolygon::nGetBoundarySize(void) const
{
   return static_cast<int>(m_VPoints.size());
}

//! Return the number of points in the up-coast profile
int CGeomCoastPolygon::nGetNumPointsUsedUpCoastProfile(void) const
{
   return m_nProfileUpCoastNumPointsUsed;
}

//! Return the number of points in the down-coast profile
int CGeomCoastPolygon::nGetNumPointsUsedDownCoastProfile(void) const
{
   return m_nProfileDownCoastNumPointsUsed;
}

//! Set the volume of seawater in the coast polygon
void CGeomCoastPolygon::SetSeawaterVolume(const double dDepth)
{
   m_dSeawaterVolume = dDepth;
}

//! Get the volume of seawater in the coast polygon
double CGeomCoastPolygon::dGetSeawaterVolume(void) const
{
   return m_dSeawaterVolume;
}

//! Adds in potential beach erosion of unconsolidated sediment (all size classes) on this polygon (m_dPotentialBeachErosionAllUncons is <= 0)
void CGeomCoastPolygon::AddPotentialErosion(double const dDepth)
{
   m_dPotentialBeachErosionAllUncons += dDepth;
}

//! Returns this timestep's total change in depth of unconsolidated sediment (all size classes) due to beach erosion on this polygon (m_dPotentialBeachErosionAllUncons <= 0)
double CGeomCoastPolygon::dGetPotentialErosion(void) const
{
   return m_dPotentialBeachErosionAllUncons;
}

//! Sets a value (must be < 0) for this timestep's erosion of fine unconsolidated sediment (beach redistribution only) on this polygon
void CGeomCoastPolygon::SetBeachErosionUnconsFine(double const dDepth)
{
   m_dBeachErosionUnconsFine = dDepth;
}

//! Returns this timestep's erosion (a value < 0) of fine unconsolidated sediment (beach redistribution only) on this polygon
double CGeomCoastPolygon::dGetBeachErosionUnconsFine(void) const
{
   return m_dBeachErosionUnconsFine;
}

//! Sets a value (must be < 0) for this timestep's erosion of sand-sized unconsolidated sediment (beach redistribution only) on this polygon. This includes sand sediment eroded during Dean profile deposition of sand sediment
void CGeomCoastPolygon::SetBeachErosionUnconsSand(double const dDepth)
{
   m_dBeachErosionUnconsSand = dDepth;
}

//! Returns this timestep's erosion (a value < 0) of sand-sized unconsolidated sediment (beach redistribution only) on this polygon. This includes sand sediment eroded during Dean profile deposition of sand sediment
double CGeomCoastPolygon::dGetBeachErosionUnconsSand(void) const
{
   return m_dBeachErosionUnconsSand;
}

//! Sets a value (must be < 0) for this timestep's erosion of coarse unconsolidated sediment (beach redistribution only) on this polygon. This includes coarse sediment eroded during Dean profile deposition of coarse sediment
void CGeomCoastPolygon::SetBeachErosionUnconsCoarse(double const dDepth)
{
   m_dBeachErosionUnconsCoarse = dDepth;
}

//! Returns this timestep's erosion (a value < 0) of coarse unconsolidated sediment (beach redistribution only) on this polygon. This includes coarse sediment eroded during Dean profile deposition of coarse sediment
double CGeomCoastPolygon::dGetBeachErosionUnconsCoarse(void) const
{
   return m_dBeachErosionUnconsCoarse;
}

//! Returns this timestep's total (all size classes) beach erosion of unconsolidated sediment on this polygon, as a -ve depth in m
double CGeomCoastPolygon::dGeBeachErosionAllUncons(void) const
{
   return m_dBeachErosionUnconsFine + m_dBeachErosionUnconsSand + m_dBeachErosionUnconsCoarse;
}

// //! Adds a depth (in m) of fine-sized unconsolidated sediment to this timestep's to-suspension movement of unconsolidated coarse sediment on this polygon
// void CGeomCoastPolygon::AddToSuspensionUnconsFine(double const dDepth)
// {
//    m_dSuspensionUnconsFine += dDepth;
// }

// //! Re-initialises this timestep's to-suspension movement of unconsolidated fine sediment on this polygon
// void CGeomCoastPolygon::SetZeroSuspensionUnconsFine(void)
// {
// m_dSuspensionUnconsFine = 0;
// }

//! Returns this timestep's to-suspension movement of fine unconsolidated sediment on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetSuspensionUnconsFine(void) const
{
   return m_dSuspensionUnconsFine;
}

//! Sets the depth of sand-sized unconsolidated sediment deposited on this polygon during this timestep (+ve, beach redistribution only)
void CGeomCoastPolygon::SetBeachDepositionUnconsSand(double const dDepth)
{
   m_dBeachDepositionUnconsSand = dDepth;
}

//! Returns the depth of sand-sized unconsolidated sediment deposited on this polygon during this timestep (+ve, beach redistribution only)
double CGeomCoastPolygon::dGetBeachDepositionUnconsSand(void) const
{
   return m_dBeachDepositionUnconsSand;
}

//! Sets the depth of coarse-sized unconsolidated sediment deposited on this polygon during this timestep (+ve, beach redistribution only)
void CGeomCoastPolygon::SetBeachDepositionUnconsCoarse(double const dDepth)
{
   m_dBeachDepositionUnconsCoarse = dDepth;
}

//! Returns the depth of coarse-sized unconsolidated sediment deposited on this polygon during this timestep (+ve, beach redistribution only)
double CGeomCoastPolygon::dGetBeachDepositionUnconsCoarse(void) const
{
   return m_dBeachDepositionUnconsCoarse;
}

//! Re-initialises this timestep's still-to-do deposition of unconsolidated sand sediment (from beach redistribution only) on this polygon
void CGeomCoastPolygon::SetZeroToDoDepositionUnconsSand(void)
{
   m_dToDoBeachDepositionUnconsSand = 0;
}

//! Adds a depth (in m) of sand-sized unconsolidated sediment to this timestep's still-to-do deposition of unconsolidated sand sediment (from beach redistribution only) on this polygon
void CGeomCoastPolygon::AddToDoBeachDepositionUnconsSand(double const dDepth)
{
   m_dToDoBeachDepositionUnconsSand += dDepth;
}

//! Returns this timestep's still-to-do deposition of sand-sized unconsolidated sediment (from beach redistribution only) on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetToDoBeachDepositionUnconsSand(void) const
{
   return m_dToDoBeachDepositionUnconsSand;
}

//! Re-initialises this timestep's still-to-do deposition of unconsolidated coarse sediment (from beach redistribution only) on this polygon
void CGeomCoastPolygon::SetZeroToDoDepositionUnconsCoarse(void)
{
   m_dToDoBeachDepositionUnconsCoarse = 0;
}

//! Adds a depth (in m) of coarse unconsolidated sediment to this timestep's still-to-do deposition of unconsolidated coarse sediment (from beach redistribution only) on this polygon
void CGeomCoastPolygon::AddToDoBeachDepositionUnconsCoarse(double const dDepth)
{
   m_dToDoBeachDepositionUnconsCoarse += dDepth;
}

//! Returns this timestep's still-to-do deposition of coarse unconsolidated sediment (from beach redistribution only) on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetToDoBeachDepositionUnconsCoarse(void) const
{
   return m_dToDoBeachDepositionUnconsCoarse;
}

//! Adds a depth (in m) of sand sediment eroded during beach Dean profile deposition
void CGeomCoastPolygon::AddBeachSandErodedDeanProfile(double const dDepth)
{
   m_dBeachSandErodedDeanProfile += dDepth;
}

//! Returns the depth (in m) of sand sediment eroded during beach Dean profile deposition
double CGeomCoastPolygon::dGetBeachSandErodedDeanProfile(void) const
{
   return m_dBeachSandErodedDeanProfile;
}

//! Adds a depth (in m) of coarse sediment eroded during beach Dean profile deposition
void CGeomCoastPolygon::AddBeachCoarseErodedDeanProfile(double const dDepth)
{
   m_dBeachCoarseErodedDeanProfile += dDepth;
}

//! Returns the depth (in m) of coarse sediment eroded during beach Dean profile deposition
double CGeomCoastPolygon::dGetBeachCoarseErodedDeanProfile(void) const
{
   return m_dBeachCoarseErodedDeanProfile;
}

//! Returns this timestep's total (all size classes) deposition and to-suspension movement of unconsolidated sediment (from beach redistribution only) on this polygon, as a +ve depth in m
double CGeomCoastPolygon::dGetBeachDepositionAndSuspensionAllUncons(void) const
{
   return m_dSuspensionUnconsFine + m_dToDoBeachDepositionUnconsSand + m_dToDoBeachDepositionUnconsCoarse;
}

//! Sets all up-coast adjacent polygons
void CGeomCoastPolygon::SetUpCoastAdjacentPolygons(vector<int> const* pnVPolygons)
{
   m_VnUpCoastAdjacentPolygon = *pnVPolygons;
}

//! Gets a single up-coast adjacent polygon
int CGeomCoastPolygon::nGetUpCoastAdjacentPolygon(int const nIndex) const
{
   // assert(nIndex < m_VnUpCoastAdjacentPolygon.size());
   return m_VnUpCoastAdjacentPolygon[nIndex];
}

//! Gets all up-coast adjacent polygons
int CGeomCoastPolygon::nGetNumUpCoastAdjacentPolygons(void) const
{
   return static_cast<int>(m_VnUpCoastAdjacentPolygon.size());
}

//! Sets all down-coast adjacent polygons
void CGeomCoastPolygon::SetDownCoastAdjacentPolygons(vector<int> const* pnVPolygons)
{
   m_VnDownCoastAdjacentPolygon = *pnVPolygons;
}

//! Gets a single down-coast adjacent polygon
int CGeomCoastPolygon::nGetDownCoastAdjacentPolygon(int const nIndex) const
{
   // assert(nIndex < m_VnDownCoastAdjacentPolygon.size());
   return m_VnDownCoastAdjacentPolygon[nIndex];
}

//! Gets all down-coast adjacent polygons
int CGeomCoastPolygon::nGetNumDownCoastAdjacentPolygons(void) const
{
   return static_cast<int>(m_VnDownCoastAdjacentPolygon.size());
}

//! Searches the list of down-coast adjacent polygons for a given polygon number
bool CGeomCoastPolygon::bDownCoastIsAlreadyPresent(int const nPoly)
{
   auto it = find(m_VnDownCoastAdjacentPolygon.begin(), m_VnDownCoastAdjacentPolygon.end(), nPoly);

   if (it == m_VnDownCoastAdjacentPolygon.end())
      // Not found
      return false;

   else
      // Found
      return true;
}

//! Searches the list of up-coast adjacent polygons for a given polygon number
bool CGeomCoastPolygon::bUpCoastIsAlreadyPresent(int const nPoly)
{
   auto it = find(m_VnUpCoastAdjacentPolygon.begin(), m_VnUpCoastAdjacentPolygon.end(), nPoly);

   if (it == m_VnUpCoastAdjacentPolygon.end())
      // Not found
      return false;

   else
      // Found
      return true;
}

//! Sets the boundary shares for all up-coast adjacent polygons
void CGeomCoastPolygon::SetUpCoastAdjacentPolygonBoundaryShares(vector<double> const* pdVShares)
{
   m_VdUpCoastAdjacentPolygonBoundaryShare = *pdVShares;
}

//! Gets the boundary shares for all up-coast adjacent polygons
double CGeomCoastPolygon::dGetUpCoastAdjacentPolygonBoundaryShare(int const nIndex) const
{
   // TODO 055 No check to see if nIndex < m_VdUpCoastAdjacentPolygonBoundaryShare.size()
   return m_VdUpCoastAdjacentPolygonBoundaryShare[nIndex];
}

//! Sets the boundary shares for all down-coast adjacent polygons
void CGeomCoastPolygon::SetDownCoastAdjacentPolygonBoundaryShares(vector<double> const* pdVShares)
{
   m_VdDownCoastAdjacentPolygonBoundaryShare = *pdVShares;
}

//! Gets the boundary shares for all down-coast adjacent polygons
double CGeomCoastPolygon::dGetDownCoastAdjacentPolygonBoundaryShare(int const nIndex) const
{
   // TODO 055 No check to see if nIndex < m_VdDownCoastAdjacentPolygonBoundaryShare.size()
   return m_VdDownCoastAdjacentPolygonBoundaryShare[nIndex];
}

//! Set the average d50 for unconsolidated sediment on this polygon
void CGeomCoastPolygon::SetAvgUnconsD50(double const dD50)
{
   m_dAvgUnconsD50 = dD50;
}

//! Get the average d50 for unconsolidated sediment on this polygon
double CGeomCoastPolygon::dGetAvgUnconsD50(void) const
{
   return m_dAvgUnconsD50;
}

//! Instantiates the pure virtual function in the abstract parent class, so that CGeomCoastPolygon is not an abstract class
void CGeomCoastPolygon::Display(void)
{
}

//! Add a circularity to this polygon
void CGeomCoastPolygon::AddCircularity(int const nPoly)
{
   m_VnCircularityWith.push_back(nPoly);
}

//! Get all circularities for this polygon
vector<int> const* CGeomCoastPolygon::VnGetCircularities(void) const
{
   return &m_VnCircularityWith;
}

//! Add to the this-iteration total of unconsolidated fine sediment from cliff collapse which goes to suspension on this polygon
void CGeomCoastPolygon::AddCliffCollapseToSuspensionFine(double const dDepth)
{
   m_dCliffCollapseTalusFineToSuspension += dDepth;
}

//! Get the this-iteration total of unconsolidated fine sediment from cliff collapse which goes to suspension on this polygon
double CGeomCoastPolygon::dGetCliffCollapseToSuspensionFine(void) const
{
   return m_dCliffCollapseTalusFineToSuspension;
}

//! Add to the this-iteration total of unconsolidated fine sediment eroded from cliff collapse on this polygon
void CGeomCoastPolygon::AddCliffCollapseErosionFine(double const dDepth)
{
   m_dCliffCollapseErosionFine += dDepth;
}

//! Get the this-iteration total of unconsolidated fine sediment eroded from cliff collapse on this polygon
double CGeomCoastPolygon::dGetCliffCollapseErosionFine(void) const
{
   return m_dCliffCollapseErosionFine;
}

//! Add to the this-iteration total of unconsolidated sand sediment from cliff collapse on this polygon, note that this does not include sand sediment eroded during Dean profile deposition of talus
void CGeomCoastPolygon::AddCliffCollapseErosionSand(double const dDepth)
{
   m_dCliffCollapseErosionSand += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment from cliff collapse on this polygon, note that this does not include sand sediment eroded during Dean profile deposition of talus
double CGeomCoastPolygon::dGetCliffCollapseErosionSand(void) const
{
   return m_dCliffCollapseErosionSand;
}

//! Add to the this-iteration total of unconsolidated coarse sediment from cliff collapse on this polygon, note that this does not include coarse sediment eroded during Dean profile deposition of talus
void CGeomCoastPolygon::AddCliffCollapseErosionCoarse(double const dDepth)
{
   m_dCliffCollapseErosionCoarse += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment from cliff collapse on this polygon, note that this does not include coarse sediment eroded during Dean profile deposition of talus
double CGeomCoastPolygon::dGetCliffCollapseErosionCoarse(void) const
{
   return m_dCliffCollapseErosionCoarse;
}

//! Add to the this-iteration total of unconsolidated sand sediment deposited from cliff collapse on this polygon
void CGeomCoastPolygon::AddCliffCollapseUnconsSandDeposition(double const dDepth)
{
   m_dCliffCollapseTalusSand += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment deposited from cliff collapse on this polygon
double CGeomCoastPolygon::dGetCliffCollapseUnconsSandDeposition(void) const
{
   return m_dCliffCollapseTalusSand;
}

//! Add to the this-iteration total of unconsolidated coarse sediment deposited from cliff collapse on this polygon
void CGeomCoastPolygon::AddCliffCollapseUnconsCoarseDeposition(double const dDepth)
{
   m_dCliffCollapseTalusCoarse += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment deposited from cliff collapse on this polygon
double CGeomCoastPolygon::dGetCliffCollapseUnconsCoarseDeposition(void) const
{
   return m_dCliffCollapseTalusCoarse;
}

//! Add to the this-iteration total of unconsolidated sand sediment eroded during deposition of cliff collapse talus as a Dean profile
void CGeomCoastPolygon::AddCliffCollapseSandErodedDeanProfile(double const dDepth)
{
   m_dCliffCollapseSandErodedDeanProfile += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment eroded during deposition of cliff collapse talus as a Dean profile
double CGeomCoastPolygon::dGetCliffCollapseSandErodedDeanProfile(void) const
{
   return m_dCliffCollapseSandErodedDeanProfile;
}

//! Add to the this-iteration total of unconsolidated coarse sediment eroded during deposition of cliff collapse talus as a Dean profile
void CGeomCoastPolygon::AddCliffCollapseCoarseErodedDeanProfile(double const dDepth)
{
   m_dCliffCollapseCoarseErodedDeanProfile += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment eroded during deposition of cliff collapse talus as a Dean profile
double CGeomCoastPolygon::dGetCliffCollapseCoarseErodedDeanProfile(void) const
{
   return m_dCliffCollapseCoarseErodedDeanProfile;
}

// //! Add to the this-iteration total of unconsolidated fine sediment moved to suspension and derived from shore platform erosion on this polygon
// void CGeomCoastPolygon::AddPlatformErosionToSuspensionUnconsFine(double const dDepth)
// {
//    m_dPlatformErosionToSuspensionFine += dDepth;
// }

// //! Get the this-iteration total of unconsolidated sand sediment moved to suspension derived from shore platform erosion on this polygon
// double CGeomCoastPolygon::dGetPlatformErosionToSuspensionUnconsFine(void) const
// {
//    return m_dPlatformErosionToSuspensionFine;
// }

//! Add to the this-iteration total of unconsolidated sand sediment derived from shore platform erosion on this polygon
void CGeomCoastPolygon::AddPlatformErosionUnconsSand(double const dDepth)
{
   m_dPlatformErosionUnconsSand += dDepth;
}

//! Get the this-iteration total of unconsolidated sand sediment derived from shore platform erosion on this polygon
double CGeomCoastPolygon::dGetPlatformErosionUnconsSand(void) const
{
   return m_dPlatformErosionUnconsSand;
}

//! Add to the this-iteration total of unconsolidated coarse sediment derived from shore platform erosion on this polygon
void CGeomCoastPolygon::AddPlatformErosionUnconsCoarse(double const dDepth)
{
   m_dPlatformErosionUnconsCoarse += dDepth;
}

//! Get the this-iteration total of unconsolidated coarse sediment derived from shore platform erosion on this polygon
double CGeomCoastPolygon::dGetPlatformErosionUnconsCoarse(void) const
{
   return m_dPlatformErosionUnconsCoarse;
}

//! Set the value of pre-existing unconsolidated fine sediment stored on this polygon
void CGeomCoastPolygon::SetPreExistingUnconsFine(double const dDepth)
{
   m_dPreExistingUnconsFine = dDepth;
}

//! Get the value of pre-existing unconsolidated fine sediment stored on this polygon
double CGeomCoastPolygon::dGetPreExistingUnconsFine(void) const
{
   return m_dPreExistingUnconsFine;
}

//! Set the value of pre-existing unconsolidated sand sediment stored on this polygon
void CGeomCoastPolygon::SetPreExistingUnconsSand(double const dDepth)
{
   m_dPreExistingUnconsSand = dDepth;
}

//! Get the value of pre-existing unconsolidated sand sediment stored on this polygon
double CGeomCoastPolygon::dGetPreExistingUnconsSand(void) const
{
   return m_dPreExistingUnconsSand;
}

//! Set the value of pre-existing unconsolidated coarse sediment stored on this polygon
void CGeomCoastPolygon::SetPreExistingUnconsCoarse(double const dDepth)
{
   m_dPreExistingUnconsCoarse = dDepth;
}

//! Get the value of pre-existing unconsolidated coarse sediment stored on this polygon
double CGeomCoastPolygon::dGetPreExistingUnconsCoarse(void) const
{
   return m_dPreExistingUnconsCoarse;
}

//! Set the value of pre-existing consolidated fine sediment stored on this polygon
void CGeomCoastPolygon::SetPreExistingConsFine(double const dDepth)
{
   m_dPreExistingConsFine = dDepth;
}

//! Get the value of pre-existing consolidated fine sediment stored on this polygon
double CGeomCoastPolygon::dGetPreExistingConsFine(void) const
{
   return m_dPreExistingConsFine;
}

//! Set the value of pre-existing consolidated sand sediment stored on this polygon
void CGeomCoastPolygon::SetPreExistingConsSand(double const dDepth)
{
   m_dPreExistingConsSand = dDepth;
}

//! Get the value of pre-existing consolidated sand sediment stored on this polygon
double CGeomCoastPolygon::dGetPreExistingConsSand(void) const
{
   return m_dPreExistingConsSand;
}

//! Set the value of pre-existing consolidated coarse sediment stored on this polygon
void CGeomCoastPolygon::SetPreExistingConsCoarse(double const dDepth)
{
   m_dPreExistingConsCoarse = dDepth;
}

//! Get the value of pre-existing consolidated coarse sediment stored on this polygon
double CGeomCoastPolygon::dGetPreExistingConsCoarse(void) const
{
   return m_dPreExistingConsCoarse;
}

//! Set the value of fine sediment on the polygon derived from sediment input events(s)
void CGeomCoastPolygon::SetSedimentInputUnconsFine(double const dDepth)
{
   m_dSedimentInputFine = dDepth;
}

//! Get the value of fine sediment on the polygon derived from sediment input events(s)
double CGeomCoastPolygon::dGetSedimentInputUnconsFine(void) const
{
   return m_dSedimentInputFine;
}

//! Set the value of sand sediment on the polygon derived from sediment input events(s)
void CGeomCoastPolygon::SetSedimentInputUnconsSand(double const dDepth)
{
   m_dSedimentInputSand = dDepth;
}

//! Get the value of sand sediment on the polygon derived from sediment input events(s)
double CGeomCoastPolygon::dGetSedimentInputUnconsSand(void) const
{
   return m_dSedimentInputSand;
}

//! Set the value of coarse sediment on the polygon derived from sediment input events(s)
void CGeomCoastPolygon::SetSedimentInputUnconsCoarse(double const dDepth)
{
   m_dSedimentInputCoarse = dDepth;
}

//! Get the value of coarse sediment on the polygon derived from sediment input events(s)
double CGeomCoastPolygon::dGetSedimentInputUnconsCoarse(void) const
{
   return m_dSedimentInputCoarse;
}

//! Appends the point cordinates (grid CRS) for a polygon vertex
void CGeomCoastPolygon::AppendVertex(CGeom2DIPoint const* pPti)
{
   m_VPtiVertices.push_back(*pPti);
}

// //! Returns the number of vertices for this polygon
// int CGeomCoastPolygon::nGetNumVertices(void) const
// {
//    return static_cast<int>(m_VPtiVertices.size());
// }

// //! Returns the point coordinates (grid CRS) for a single vertex of this polygon
// CGeom2DIPoint CGeomCoastPolygon::PtiGetVertex(int const nIndex) const
// {
//    // Note no check to see if nUndex < m_VPtiVertices.size()
//    return m_VPtiVertices[nIndex];
// }

CGeom2DIPoint CGeomCoastPolygon::PtiGetFillStartPoint(void)
{
   int const nVertices = static_cast<int>(m_VPtiVertices.size());
   double dXTot = 0;
   double dYTot = 0;

   for (int n = 0; n < nVertices; n++)
   {
      dXTot += m_VPtiVertices[n].nGetX();
      dYTot += m_VPtiVertices[n].nGetY();
   }

   return CGeom2DIPoint(nRound(dXTot / nVertices), nRound(dYTot / nVertices));
}
