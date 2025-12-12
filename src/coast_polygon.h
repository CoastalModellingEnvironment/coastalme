/*!
   \class CGeomCoastPolygon
   \brief Geometry class used for coast polygon objects
   \details TODO 001 This is a more detailed description of the CRWCoast class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file coast_polygon.h
   \brief Contains CGeomCoastPolygon definitions

*/

#ifndef COASTPOLYGON_H
#define COASTPOLYGON_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "2d_point.h"
#include "2di_point.h"
#include "2d_shape.h"

class CGeomCoastPolygon : public CA2DShape
{
 private:
   //! Is the movement of unconsolidated sediment on this polygon down-coast (i.e. along the coast in the direction of increasing coastline point numbers) during this iteration?
   bool m_bUnconsSedimentMovementDownCoastThisIter;

   // Does the polygon meet at a point at its seaward end? (is it roughly triangular?)
   // bool m_bIsPointedSeaward;

   //! Is this polygon at the end of the coastline?
   bool m_bCoastEndPolygon;

   //! Is this polygon at the start of the coastline?
   bool m_bCoastStartPolygon;

   //! This-coast-only number of this polygon
   int m_nCoastID;

   //! The point on this polygon's coastline segment with maximum concave curvature, roughly at the middle of the coastline segment
   int m_nCoastNode;

   //! The normal profile which bounds the polygon in the up-coast direction
   int m_nProfileUpCoast;

   //! The normal profile which bounds the polygon in the down-coast direction
   int m_nProfileDownCoast;

   //! The number of points from the up-coast normal which are part of this polygon (less than the normal's full length if the polygon is triangular)
   int m_nProfileUpCoastNumPointsUsed;

   //! The number of points from the down-coast normal which are part of this polygon (less than the normal's full length if the polygon is triangular)
   int m_nProfileDownCoastNumPointsUsed;

   //! The number of cells in the polygon
   int m_nNumCells;

   //! The average d50 of unconsolidated sediment on this polygon
   double m_dAvgUnconsD50;

   //! The volume (m3) of sea water within the polygon
   double m_dSeawaterVolume;

   // Note: all sediment depths are in m, and here cover the area of a single raster cell: to convert to a volume, multiply by m_dCellArea

   //! Potential (ignoring supply-limitation) erosion of unconsolidated sediment (all size classes) as a depth during this timestep (-ve), beach redistribution only
   double m_dPotentialBeachErosionAllUncons;

   //! Erosion (considering supply-limitation) of fine-sized unconsolidated sediment as a depth this timestep (-ve), beach redistribution only
   double m_dBeachErosionUnconsFine;

   //! Erosion (considering supply-limitation) of sand-sized unconsolidated sediment as a depth this timestep (-ve), beach redistribution only. This includes sand sediment eroded during Dean profile deposition of sand sediment
   double m_dBeachErosionUnconsSand;

   //! Erosion (considering supply-limitation) of coarse-sized unconsolidated sediment as a depth this timestep (-ve), beach redistribution only. This includes coarse sediment eroded during Dean profile deposition of coarse sediment
   double m_dBeachErosionUnconsCoarse;

   //! Depth of sand-sized unconsolidated sediment deposited on this polygon during this timestep (+ve, beach redistribution only)
   double m_dBeachDepositionUnconsSand;

   //! Depth of coarse-sized unconsolidated sediment deposited on this polygon during this timestep (+ve, beach redistribution only)
   double m_dBeachDepositionUnconsCoarse;

   //! To-suspension movement of fine-sized sediment as a depth this timestep (+ve) ALL PROCESSES TODO *** CHECK *** REPLACE, INSTEAD USE TOTALS FROM TO-SUSPENSION FINE FROM (a) PLATFORM EROSION (b) CLIFF COLLAPSE (c) BEACH EROSION
   double m_dSuspensionUnconsFine;

   //! Still-to-do depth (m) of sand-sized unconsolidated sediment to be deposited this timestep (+ve), beach redistribution only
   double m_dToDoBeachDepositionUnconsSand;

   //! Still-to-do depth (m) of coarse-sized unconsolidated sediment to be deposited this timestep (+ve), beach redistribution only
   double m_dToDoBeachDepositionUnconsCoarse;

   //! Depth of sand unconsolidated sediment eroded during beach deposition as a Dean profile
   double m_dBeachSandErodedDeanProfile;

   //! Depth of coarse unconsolidated sediment eroded during beach deposition as a Dean profile
   double m_dBeachCoarseErodedDeanProfile;

   //! Depth of eroded fine sediment from cliff collapse (is always equal to m_dCliffCollapseToSuspensionFine)
   double m_dCliffCollapseErosionFine;

   //! Depth of eroded sand sediment from cliff collapse, note that this does not include sand sediment eroded during Dean profile deposition of talus
   double m_dCliffCollapseErosionSand;

   //! Depth of eroded coarse sediment from cliff collapse, note that this does not include coarse sediment eroded during Dean profile deposition of talus
   double m_dCliffCollapseErosionCoarse;

   //! Depth of unconsolidated fine sediment which goes to suspension from cliff collapse
   double m_dCliffCollapseTalusFineToSuspension;

   //! Depth of unconsolidated sand talus from cliff collapse
   double m_dCliffCollapseTalusSand;

   //! Depth of unconsolidated coarse talus from cliff collapse
   double m_dCliffCollapseTalusCoarse;

   //! Depth of sand unconsolidated sediment eroded during deposition of cliff collapse talus as a Dean profile
   double m_dCliffCollapseSandErodedDeanProfile;

   //! Depth of coarse unconsolidated sediment eroded during deposition of cliff collapse talus as a Dean profile
   double m_dCliffCollapseCoarseErodedDeanProfile;

   //! Depth of fine sediment moved to suspension from shore platform erosion
   double m_dPlatformErosionToSuspensionFine;

   //! Depth of unconsolidated sand sediment from shore platform erosion
   double m_dPlatformErosionUnconsSand;

   //! Depth of unconsolidated coarse sediment from shore platform erosion
   double m_dPlatformErosionUnconsCoarse;

   //! Depth of pre-existing unconsolidated fine sediment
   double m_dPreExistingUnconsFine;

   //! Depth of pre-existing unconsolidated sand sediment
   double m_dPreExistingUnconsSand;

   //! Depth of pre-existing unconsolidated coarse sedimet
   double m_dPreExistingUnconsCoarse;

   //! Depth of pre-existing consolidated fine sediment
   double m_dPreExistingConsFine;

   //! Depth of pre-existing consolidated sand sediment
   double m_dPreExistingConsSand;

   //! Depth of pre-existing consolidated coarse sediment
   double m_dPreExistingConsCoarse;

   //! Depth of fine sediment added from this-iteration sediment input event(s)
   double m_dSedimentInputFine;

   //! Depth of sand sediment added from this-iteration sediment input event(s)
   double m_dSedimentInputSand;

   //! Depth of coarse sediment added from this-iteration sediment input event(s)
   double m_dSedimentInputCoarse;

   //! Coast polygon length
   double m_dLength;

   //! Coordinates of the coast node cell (raster grid CRS)
   CGeom2DIPoint m_PtiNode;

   //! Coordinates of the cell (raster grid CRS) which is at other (seaward) end of the polygon
   CGeom2DIPoint m_PtiAntinode;

   //! The ID(s) of the up-coast adjacent polygon(s)
   vector<int> m_VnUpCoastAdjacentPolygon;

   //! The ID(s) of the down-coast adjacent polygon(s)
   vector<int> m_VnDownCoastAdjacentPolygon;

   //! If this polygon has a circular unconsolidated-sediment-movement relationship with one or more other polygons, the coast-only numbers of these polygons
   vector<int> m_VnCircularityWith;

   //! The boundary share(s) (0 to 1) with adjacent up-coast polygon(s)
   vector<double> m_VdUpCoastAdjacentPolygonBoundaryShare;

   //! The boundary share(s) (0 to 1) with adjacent up-coast polygon(s)
   vector<double> m_VdDownCoastAdjacentPolygonBoundaryShare;

   //! The polygon's vertices (not all of them, just the ends of the profile sides), used to calculate the polygon's centroid for filling it
   vector<CGeom2DIPoint> m_VPtiVertices;

 protected:
 public:
   CGeomCoastPolygon(int const, int const, int const, int const, vector<CGeom2DPoint> const*, int const, int const, CGeom2DIPoint const*, CGeom2DIPoint const*, bool const, bool const);
   ~CGeomCoastPolygon(void) override;

   void SetDownCoastThisIter(bool const);
   bool bDownCoastThisIter(void) const;

   // void SetCoastEndPolygon(void);
   bool bIsCoastEndPolygon(void) const;
   // void SetCoastStartPolygon(void);
   bool bIsCoastStartPolygon(void) const;

   int nGetPolygonCoastID(void) const;

   // void SetCoastNode(int const);
   int nGetNodeCoastPoint(void) const;
   CGeom2DIPoint* pPtiGetNode(void);
   CGeom2DIPoint* pPtiGetAntiNode(void);

   void SetLength(double const);
   double dGetLength(void) const;

   // void SetNotPointed(void);
   // bool bIsPointed(void) const;

   void SetNumCellsInPolygon(int const);
   // int nGetNumCellsinPolygon(void) const;

   int nGetUpCoastProfile(void) const;
   int nGetDownCoastProfile(void) const;

   // void SetBoundary(vector<CGeom2DPoint> const*);
   // vector<CGeom2DPoint>* pPtVGetBoundary(void);
   CGeom2DPoint* pPtGetBoundaryPoint(int const);
   int nGetBoundarySize(void) const;

   int nGetNumPointsUsedUpCoastProfile(void) const;
   int nGetNumPointsUsedDownCoastProfile(void) const;

   void SetSeawaterVolume(const double);
   double dGetSeawaterVolume(void) const;

   void AddPotentialErosion(double const);
   double dGetPotentialErosion(void) const;

   void SetBeachErosionUnconsFine(double const);
   double dGetBeachErosionUnconsFine(void) const;
   void SetBeachErosionUnconsSand(double const);
   double dGetBeachErosionUnconsSand(void) const;
   void SetBeachErosionUnconsCoarse(double const);
   double dGetBeachErosionUnconsCoarse(void) const;
   double dGeBeachErosionAllUncons(void) const;

   void SetBeachDepositionUnconsSand(double const);
   double dGetBeachDepositionUnconsSand(void) const;
   void SetBeachDepositionUnconsCoarse(double const);
   double dGetBeachDepositionUnconsCoarse(void) const;

   // void AddToSuspensionUnconsFine(double const);
   // void SetZeroSuspensionUnconsFine(void);
   double dGetSuspensionUnconsFine(void) const;

   void SetZeroToDoDepositionUnconsSand(void);
   void AddToDoBeachDepositionUnconsSand(double const);
   double dGetToDoBeachDepositionUnconsSand(void) const;
   void SetZeroToDoDepositionUnconsCoarse(void);
   void AddToDoBeachDepositionUnconsCoarse(double const);
   double dGetToDoBeachDepositionUnconsCoarse(void) const;
   double dGetBeachDepositionAndSuspensionAllUncons(void) const;

   void AddBeachSandErodedDeanProfile(double const);
   double dGetBeachSandErodedDeanProfile(void) const;
   void AddBeachCoarseErodedDeanProfile(double const);
   double dGetBeachCoarseErodedDeanProfile(void) const;

   void SetUpCoastAdjacentPolygons(vector<int> const*);
   int nGetUpCoastAdjacentPolygon(int const) const;
   int nGetNumUpCoastAdjacentPolygons(void) const;

   void SetDownCoastAdjacentPolygons(vector<int> const*);
   int nGetDownCoastAdjacentPolygon(int const) const;
   int nGetNumDownCoastAdjacentPolygons(void) const;

   void SetUpCoastAdjacentPolygonBoundaryShares(vector<double> const*);
   double dGetUpCoastAdjacentPolygonBoundaryShare(int const) const;

   void SetDownCoastAdjacentPolygonBoundaryShares(vector<double> const*);
   double dGetDownCoastAdjacentPolygonBoundaryShare(int const) const;

   bool bDownCoastIsAlreadyPresent(int const);
   bool bUpCoastIsAlreadyPresent(int const);

   // int nGetPointInPolygonSearchStartPoint(void) const;

   void SetAvgUnconsD50(double const);
   double dGetAvgUnconsD50(void) const;

   void Display(void) override;

   void AddCircularity(int const);
   vector<int> const* VnGetCircularities(void) const;

   void AddCliffCollapseErosionFine(double const);
   double dGetCliffCollapseErosionFine(void) const;
   void AddCliffCollapseErosionSand(double const);
   double dGetCliffCollapseErosionSand(void) const;
   void AddCliffCollapseErosionCoarse(double const);
   double dGetCliffCollapseErosionCoarse(void) const;

   void AddCliffCollapseToSuspensionFine(double const);
   double dGetCliffCollapseToSuspensionFine(void) const;
   void AddCliffCollapseUnconsSandDeposition(double const);
   double dGetCliffCollapseUnconsSandDeposition(void) const;
   void AddCliffCollapseUnconsCoarseDeposition(double const);
   double dGetCliffCollapseUnconsCoarseDeposition(void) const;

   void AddCliffCollapseFineErodedDeanProfile(double const);
   double dGetCliffCollapseFineErodedDeanProfile(void) const;
   void AddCliffCollapseSandErodedDeanProfile(double const);
   double dGetCliffCollapseSandErodedDeanProfile(void) const;
   void AddCliffCollapseCoarseErodedDeanProfile(double const);
   double dGetCliffCollapseCoarseErodedDeanProfile(void) const;

   // void AddPlatformErosionToSuspensionUnconsFine(double const);
   // double dGetPlatformErosionToSuspensionUnconsFine(void) const;
   void AddPlatformErosionUnconsSand(double const);
   double dGetPlatformErosionUnconsSand(void) const;
   void AddPlatformErosionUnconsCoarse(double const);
   double dGetPlatformErosionUnconsCoarse(void) const;

   void SetPreExistingUnconsFine(double const);
   double dGetPreExistingUnconsFine(void) const;
   void SetPreExistingUnconsSand(double const);
   double dGetPreExistingUnconsSand(void) const;
   void SetPreExistingUnconsCoarse(double const);
   double dGetPreExistingUnconsCoarse(void) const;

   void SetPreExistingConsFine(double const);
   double dGetPreExistingConsFine(void) const;
   void SetPreExistingConsSand(double const);
   double dGetPreExistingConsSand(void) const;
   void SetPreExistingConsCoarse(double const);
   double dGetPreExistingConsCoarse(void) const;

   void SetSedimentInputUnconsFine(double const);
   double dGetSedimentInputUnconsFine(void) const;
   void SetSedimentInputUnconsSand(double const);
   double dGetSedimentInputUnconsSand(void) const;
   void SetSedimentInputUnconsCoarse(double const);
   double dGetSedimentInputUnconsCoarse(void) const;

   void AppendVertex(CGeom2DIPoint const*);
   // int nGetNumVertices(void) const;
   // CGeom2DIPoint PtiGetVertex(int const) const;

   CGeom2DIPoint PtiGetFillStartPoint(void);
};
#endif // COASTPOLYGON_H
