/*!
   \class CGeomProfile
   \brief Geometry class used to represent coast profile objects
   \details TODO 001 This is a more detailed description of the CGeomProfile class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file profile.h
   \brief Contains CGeomProfile definitions
*/

#ifndef PROFILE_H
#define PROFILE_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <vector>
using std::vector;

#include "2d_point.h"
#include "2di_point.h"
#include "multi_line.h"
#include "raster_grid.h"

class CGeomProfile : public CGeomMultiLine
{
 private:
   //! Is this a start-of-coast profile?
   bool m_bStartOfCoast;

   //! Is this an end-of-coast profile?
   bool m_bEndOfCoast;

   //! Has this profile encountered a CShore problem?
   bool m_bCShoreProblem;

   //! Has this profile hit land?
   bool m_bHitLand;

   //! Has this profile hit an intervention?
   bool m_bHitIntervention;

   //! Has this profile hit a coastline?
   bool m_bHitCoast;

   //! Is this profile too short?
   bool m_bTooShort;

   //! Has this profile been truncated by hitting another profile from the same coast?
   bool m_bTruncatedSameCoast;

   //! Has this profile been truncated by hitting another profile from a different coast?
   bool m_bTruncatedDifferentCoast;

   //! Has this profile hit another profile?
   bool m_bHitAnotherProfile;

   //! Is this an intervention profile?
   bool m_bIntervention;

   //! The coast from which this profile projects
   int m_nCoast;

   //! The coastline point at which this profile hits the coast (not necessarily coincident wih the profile start cell)
   int m_nCoastPoint;

   //! The this-coast ID of the profile (note that a profile in a different coast may have the same ID as this profile)
   int m_nProfileID;

   //! The wave height at the end of the profile
   double m_dDeepWaterWaveHeight;

   //! The wave orientation at the end of the profile
   double m_dDeepWaterWaveAngle;

   //! The wave period at the end of the profile
   double m_dDeepWaterWavePeriod;

   //! Pointer to the adjacent up-coast profile (may be an invalid profile)
   CGeomProfile* m_pUpCoastAdjacentProfile;

   //! Pointer to the adjacent down-coast profile (may be an invalid profile)
   CGeomProfile* m_pDownCoastAdjacentProfile;

   //! In the grid CRS, the integer coordinates of the cells 'under' this profile, point zero is the same as 'cell marked as coastline' in coast object
   vector<CGeom2DIPoint> m_VCellInProfile;

   //! In external CRS, the coords of cells 'under' this profile (has the same length as m_VCellInProfile)
   // vector<CGeom2DPoint> m_VCellInProfileExtCRS;

   // Is this profile point part of a multi-line? (Has have the same length as m_VCellInProfile)
   // vector<bool> m_bVShared;

 protected:
 public:
   explicit CGeomProfile(int const, int const, int const, bool const);
   ~CGeomProfile(void) override;

   int nGetCoastID(void) const;
   int nGetProfileID(void) const;
   int nGetCoastPoint(void) const;

   CGeom2DIPoint* pPtiGetStartPoint(void);
   CGeom2DIPoint* pPtiGetEndPoint(void);

   void SetStartOfCoast(bool const);
   bool bStartOfCoast(void) const;
   void SetEndOfCoast(bool const);
   bool bEndOfCoast(void) const;
   bool bIsGridEdge(void) const;

   void SetCShoreProblem(bool const);
   bool bCShoreProblem(void) const;

   void SetHitLand(bool const);
   bool bHitLand(void) const;
   void SetHitIntervention(bool const);
   bool bHitIntervention(void) const;
   void SetHitCoast(bool const);
   bool bHitCoast(void) const;
   void SetTooShort(bool const);
   bool bTooShort(void) const;
   void SetTruncatedSameCoast(bool const);
   bool bTruncatedSameCoast(void) const;
   void SetTruncatedDifferentCoast(bool const);
   bool bTruncatedDifferentCoast(void) const;
   void SetHitAnotherProfile(bool const);
   bool bHitAnotherProfile(void) const;

   bool bProfileOK(void) const;
   bool bProfileOKIncTruncated(void) const;
   bool bOKIncStartAndEndOfCoast(void) const;
   // bool bOKIncStartOfCoast(void) const;

   void SetPointsInProfile(vector<CGeom2DPoint> const*);
   void SetPointInProfile(int const, double const, double const);
   void AppendPointInProfile(double const, double const);
   void AppendPointInProfile(CGeom2DPoint const*);
   void TruncateProfile(int const);
   // void TruncateAndSetPointInProfile(int const, double const, double const);
   bool bInsertIntersection(double const, double const, int const);
   // void ShowProfile(void) const;
   int nGetProfileSize(void) const;
   CGeom2DPoint* pPtGetPointInProfile(int const);
   vector<CGeom2DPoint> PtVGetThisPointAndAllAfter(int const);
   // void RemoveLineSegment(int const);
   bool bIsPointInProfile(double const, double const);
   bool bIsPointInProfile(double const, double const, int&);
   // int nFindInsertionLineSeg(double const, double const);

   // void AppendPointShared(bool const);
   // bool bPointShared(int const) const;

   void SetUpCoastAdjacentProfile(CGeomProfile*);
   // CGeomProfile* pGetUpCoastAdjacentProfile(void) const;
   void SetDownCoastAdjacentProfile(CGeomProfile*);
   CGeomProfile* pGetDownCoastAdjacentProfile(void) const;

   void AppendCellInProfile(CGeom2DIPoint const*);
   void AppendCellInProfile(int const, int const);
   void SetCellsInProfile(vector<CGeom2DIPoint> const*);
   vector<CGeom2DIPoint>* pPtiVGetCellsInProfile(void);
   CGeom2DIPoint* pPtiGetCellInProfile(int const);
   int nGetNumCellsInProfile(void) const;
   // int nGetIndexOfCellInProfile(int const, int const);
   CGeom2DIPoint* pPtiGetLastCellInProfile(void);
   CGeom2DIPoint* pPtiGetFirstCellInProfile(void);

   int nGetCellGivenDepth(CGeomRasterGrid const*, double const);

   void SetProfileDeepWaterWaveHeight(double const);
   double dGetProfileDeepWaterWaveHeight(void) const;

   void SetProfileDeepWaterWaveAngle(double const);
   double dGetProfileDeepWaterWaveAngle(void) const;

   void SetProfileDeepWaterWavePeriod(double const);
   double dGetProfileDeepWaterWavePeriod(void) const;

   bool bIsIntervention(void) const;
};
#endif // PROFILE_H
