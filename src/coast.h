/*!
   \class CRWCoast
   \brief Real-world class used to represent coastline objects
   \details TODO 001 This is a more detailed description of the CRWCoast class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file coast.h
   \brief Contains CRWCoast definitions
*/

#ifndef COAST_H
#define COAST_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include "simulation.h"
#include "profile.h"
#include "cell.h"
#include "coast_landform.h"
#include "coast_polygon.h"
#include "line.h"
#include "i_line.h"
#include "2d_point.h"
#include "2di_point.h"

class CGeomProfile;      // Forward declarations
class CACoastLandform;
class CGeomCoastPolygon;

class CRWCoast
{
   private:
   //! Direction of the sea from the coastline, travelling down-coast (i.e. in direction of increasing coast point indices)
   int m_nSeaHandedness;

   //! The edge from which the coast starts
   int m_nStartEdge;

   //! The edge at which the coast ends
   int m_nEndEdge;

   //! The mean of the coast's detailed curvature
   double m_dCurvatureDetailedMean;

   //! The standard deviation of the coast's detailed curvature
   double m_dCurvatureDetailedSTD;

   //! The mean of the coast's smoothed curvature
   double m_dCurvatureSmoothMean;

   //! The standard deviaton of the coast's smoothed curvature
   double m_dCurvatureSmoothSTD;

   //! A pointer to the CSimulation object
   CSimulation* m_pSim;

   //! Smoothed line of points (external CRS) giving the plan view of the vector coast
   CGeomLine m_LCoastlineExtCRS;

   //! Line of points (external CRS) giving the plan view of the vector flood of wave setup
   CGeomLine m_LFloodWaveSetupExtCRS;

   //! Line of points (external CRS) giving the plan view of the vector flood of wave setup + surge
   CGeomLine m_LFloodWaveSetupSurgeExtCRS;

   //! Line of points (external CRS) giving the plan view of the vector flood of wave setup + surge + runup
   CGeomLine m_LFloodWaveSetupSurgeRunupExtCRS;

   //! Unsmoothed integer x-y coordinates (grid CRS) of the cell marked as coastline for each point on the vector coastline. Note that where there is a coast-normal profile, this is the same as point zero in the profile coordinates
   CGeomILine m_ILCellsMarkedAsCoastline;

   //! Distance of breaking (in cells), at each point on m_LCoastlineExtCRS
   vector<int> m_VnBreakingDistance;

   //! At every point on m_LCoastlineExtCRS: INT_NODATA if no nodepoint there, otherwise the node (point of greatest concave curvature) number for a coast polygon
   vector<int> m_VnPolygonNode;

   //! Detailed curvature at each point on m_LCoastlineExtCRS
   vector<double> m_VdCurvatureDetailed;

   //! Smoothed curvature at each point on m_LCoastlineExtCRS
   vector<double> m_VdCurvatureSmooth;

   //! The deep water wave height at the end of a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdDeepWaterWaveHeight;

   //! The deep water wave orientation at the end of a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdDeepWaterWaveAngle;

   //! The deep water wave period at the end of a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdDeepWaterWavePeriod;

   //! The breaking wave height on a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdBreakingWaveHeight;

   //! The wave setup on a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdWaveSetupSurge;

   // The storm surge on a normal drawn from each point on m_LCoastlineExtCRS
   // vector<double> m_VdStormSurge;

   //! The run-up on a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdRunUp;

   //! The wave height at coast point on a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdCoastWaveHeight;

   //! The breaking wave orientation on a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdBreakingWaveAngle;

   //! The depth of breaking on a normal drawn from each point on m_LCoastlineExtCRS
   vector<double> m_VdDepthOfBreaking;

   //! As in the COVE model, this is the orientation alongshore energy/sediment movement; a +ve flux is in direction of increasing indices along coast. At each point on m_LCoastlineExtCRS
   vector<double> m_VdFluxOrientation;

   //! Wave energy at each point on m_LCoastlineExtCRS
   vector<double> m_VdWaveEnergyAtBreaking;

   //! Pointer to a coast landform object, at each point on the coastline
   vector<CACoastLandform*> m_pVLandform;

   //! Pointers to coast-normal profile objects, one for each point on the coastline (is null for most coastline points)
   vector<CGeomProfile*> m_pVNormalProfileDownAllCoastpointSeq;

   // These do not have the same length as m_LCoastlineExtCRS

   //! Coast-normal profile objects, in sequence of creation (which is the same as nGetProfileID() sequence)
   vector<CGeomProfile*> m_pVProfile;

   //! Pointers to coastline-normal profile objects, in along-coastline sequence
   vector<CGeomProfile*> m_pVProfileDownCoastSeq;

   //! This coast's polygon objects, in down-coast sequence
   vector<CGeomCoastPolygon*> m_pVCoastPolygon;

   //! Lines which comprise the edge of a shadow zone, ext CRS
   vector<CGeomLine> m_LShadowBoundary;

   //! Lines which comprise the edge of a downdrift zone, ext CRS
   vector<CGeomLine> m_LShadowDowndriftBoundary;

   protected:
   public:
   explicit CRWCoast(CSimulation*);
   ~CRWCoast(void);

   // CSimulation* pGetSim(void) const;

   void SetSeaHandedness(int const);
   int nGetSeaHandedness(void) const;

   void SetStartEdge(int const);
   int nGetStartEdge(void) const;

   void SetEndEdge(int const);
   int nGetEndEdge(void) const;

   void SetCoastlineExtCRS(CGeomLine const*);
   CGeomLine* pLGetCoastlineExtCRS(void);
   // CGeomLine* pLGetFloodWaveSetupExtCRS(void);
   // void SetFloodWaveSetupPointExtCRS(CGeomLine const*);
   // void SetFloodWaveSetupSurgePointExtCRS(CGeomLine const*);
   // void SetFloodWaveSetupSurgeRunupPointExtCRS(CGeomLine const*);
   CGeom2DPoint* pPtGetCoastlinePointExtCRS(int const);
   // CGeom2DPoint* pPtGetFloodWaveSetupPointExtCRS(int const);
   // CGeom2DPoint* pPtGetFloodWaveSetupSurgePointExtCRS(int const);
   // CGeom2DPoint* pPtGetFloodWaveSetupSurgeRunupPointExtCRS(int const);

   int nGetCoastlineSize(void) const;
   // double dGetCoastlineSegmentLength(int const, int const);
   // double dGetCoastlineLengthSoFar(int const);
   // void DisplayCoastline(void);

   void SetCoastlineGridCRS(CGeomILine const*);
   CGeomILine* pILGetCoastlineGridCRS(void);
   // void AppendCellMarkedAsCoastline(CGeom2DIPoint const*);
   // void AppendCellMarkedAsCoastline(int const, int const);
   CGeom2DIPoint* pPtiGetCellMarkedAsCoastline(int const);
   // int nGetNCellsMarkedAsCoastline(void) const;
   int nGetCoastPointGivenCell(CGeom2DIPoint*);

   double dGetDetailedCurvature(int const) const;
   void SetDetailedCurvature(int const, double const);
   vector<double>* pVGetDetailedCurvature(void);
   double dGetSmoothCurvature(int const) const;
   void SetSmoothCurvature(int const, double const);
   vector<double>* pVGetSmoothCurvature(void);
   void SetDetailedCurvatureMean(double const);
   // double dGetDetailedCurvatureMean(void) const;
   void SetDetailedCurvatureSTD(double const);
   // double dGetDetailedCurvatureSTD(void) const;
   void SetSmoothCurvatureMean(double const);
   // double dGetSmoothCurvatureMean(void) const;
   void SetSmoothCurvatureSTD(double const);
   // double dGetSmoothCurvatureSTD(void) const;

   void AppendProfile(CGeomProfile*);
   CGeomProfile* pGetProfile(int const);
   // CGeomProfile* pGetLastProfile(void);
   // void ReplaceProfile(int const, vector<CGeom2DPoint> const*);
   int nGetNumProfiles(void) const;
   void CreateProfileDownCoastIndex(void);
   void InsertProfilesInProfileCoastPointIndex(void);

   // CGeomProfile* pGetDownCoastProfile(CGeomProfile const* pProfile);
   // CGeomProfile* pGetDownCoastProfileNotIncLastProfile(CGeomProfile const* pProfile);
   // CGeomProfile* pGetUpCoastProfile(CGeomProfile const* pProfile);

   void CreateProfilesAtCoastPoints(void);
   void SetProfileAtCoastPoint(int const, CGeomProfile* const);
   bool bIsProfileAtCoastPoint(int const) const;
   CGeomProfile* pGetProfileAtCoastPoint(int const) const;
   CGeomProfile* pGetProfileWithDownCoastSeq(int const) const;
   CGeomProfile* pGetProfileWithUpCoastSeq(int const) const;

   CGeomCoastPolygon* pPolyCreateAndAppendPolygon(int const, int const, CGeom2DIPoint const*, CGeom2DIPoint const*, int const, int const, vector<CGeom2DPoint> const*, int const, int const, bool const, bool const);
   int nGetNumPolygons(void) const;
   CGeomCoastPolygon* pGetPolygon(int const);

   void SetPolygonNode(int const, int const);
   int nGetPolygonNode(int const) const;

   void SetCoastDeepWaterWaveHeight(int const, double const);
   // double dGetCoastDeepWaterWaveHeight(int const) const;

   void SetCoastDeepWaterWaveAngle(int const, double const);
   double dGetCoastDeepWaterWaveAngle(int const) const;

   void SetCoastDeepWaterWavePeriod(int const, double const);
   double dGetCoastDeepWaterWavePeriod(int const) const;

   void SetBreakingWaveHeight(int const, double const);
   double dGetBreakingWaveHeight(int const) const;

   void SetCoastWaveHeight(int const, double const);
   double dGetCoastWaveHeight(int const) const;

   void SetBreakingWaveAngle(int const, double const);
   double dGetBreakingWaveAngle(int const) const;

   void SetWaveSetupSurge(int const, double const);
   double dGetWaveSetupSurge(int const) const;

   // void SetStormSurge(int const, double const);
   // double dGetStormSurge(int const) const;

   void SetRunUp(int const, double const);
   double dGetRunUp(int const) const;

   double dGetLevel(int const, int const) const;

   void SetDepthOfBreaking(int const, double const);
   double dGetDepthOfBreaking(int const) const;

   void SetBreakingDistance(int const, int const);
   int nGetBreakingDistance(int const) const;

   void SetFluxOrientation(int const, double const);
   double dGetFluxOrientation(int const) const;

   void SetWaveEnergyAtBreaking(int const, double const);
   double dGetWaveEnergyAtBreaking(int const) const;

   void AppendCoastLandform(CACoastLandform*);
   CACoastLandform* pGetCoastLandform(int const);

   // void AppendPolygonLength(const double);
   // double dGetPolygonLength(int const) const;

   int nGetNumShadowBoundaries(void) const;
   void AppendShadowBoundary(CGeomLine const*);
   CGeomLine* pGetShadowBoundary(int const);

   int nGetNumShadowDowndriftBoundaries(void) const;
   void AppendShadowDowndriftBoundary(CGeomLine const*);
   CGeomLine* pGetShadowDowndriftBoundary(int const);
};
#endif      // COAST_H
