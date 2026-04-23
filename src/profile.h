/*!
   \class CGeomProfile
   \brief Geometry class used to represent coast profile objects
   \details TODO 001 This is a more detailed description of the CGeomProfile class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
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
   bool m_bIsStartOfCoast;

   //! Is this an end-of-coast profile?
   bool m_bIsEndOfCoast;

   //! Is this an intervention profile?
   bool m_bIntervention;

   //! Does this profile have a CShore problem?
   bool m_bCShoreProblem;

   //! The coast from which this profile projects
   int m_nCoast;

   //! The coastline point at which this profile hits the coast (not necessarily coincident wih the profile start cell)
   int m_nCoastPoint;

   //! The this-coast ID of the profile (note that a profile belonging to a different coast may have the same ID as this profile)
   int m_nProfileID;

   //! The profile's status
   int m_nProfileStatus;

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

   //! In the grid CRS, integer coordinates of the cells 'under' this profile. Point zero is the same as 'cell marked as coastline' in coast object. Vector is not filled if the profile is invalid
   vector<CGeom2DIPoint> m_VCellInProfile;

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
   bool bIsStartOfCoast(void) const;
   void SetEndOfCoast(bool const);
   bool bIsEndOfCoast(void) const;
   bool bIsStartOrEndOfCoast(void) const;
   bool bIsIntervention(void) const;
   void SetCShoreProblem(void);
   bool bHasCShoreProblem(void) const;

   void SetProfileStatus(int const);
   int nGetProfileStatus(void) const;
   bool bProfileOK(void) const;

   void SetPointsInProfile(vector<CGeom2DPoint> const*);
   bool bSetPointInProfile(int const, double const, double const);
   void AppendPointInProfile(double const, double const);
   void AppendPointInProfile(CGeom2DPoint const*);
   void TruncateProfile(int const);
   // void TruncateAndbSetPointInProfile(int const, double const, double const);
   bool bInsertIntersection(double const, double const, int const);
   // void ShowProfile(void) const;
   int nGetProfileSize(void) const;
   CGeom2DPoint* pPtGetPointInProfile(int const);
   vector<CGeom2DPoint> PtVGetThisPointAndAllAfter(int const);
   // void RemoveLineSegment(int const);
   bool bIsPointInProfile(double const, double const);
   bool bIsPointInProfile(double const, double const, int&);
   // int nFindInsertionLineSeg(double const, double const);

   void SetUpCoastAdjacentProfile(CGeomProfile*);
   CGeomProfile* pGetUpCoastAdjacentProfile(void) const;
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
};
#endif // PROFILE_H
