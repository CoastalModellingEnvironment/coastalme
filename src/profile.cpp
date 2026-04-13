/*!
   \file profile.cpp
   \brief CGeomProfile routines
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
#include <assert.h>

#include <cstdio>

#include <vector>
using std::vector;

#include <algorithm>
using std::find;

#include "cme.h"
#include "cell.h"
#include "2d_point.h"
#include "2di_point.h"
#include "multi_line.h"
#include "raster_grid.h"
#include "profile.h"

//! Constructor with initialisation list
CGeomProfile::CGeomProfile(int const nCoast, int const nCoastPoint, int const nProfileID, bool const bIntervention)
    : m_bIsStartOfCoast(false),
      m_bIsEndOfCoast(false),
      m_bIntervention(bIntervention),
      m_bCShoreProblem(false),
      m_nCoast(nCoast),
      m_nCoastPoint(nCoastPoint),
      m_nProfileID(nProfileID),
      m_nProfileStatus(PROFILE_STATUS_OK),
      m_dDeepWaterWaveHeight(0),
      m_dDeepWaterWaveAngle(0),
      m_dDeepWaterWavePeriod(0),
      m_pUpCoastAdjacentProfile(NULL),
      m_pDownCoastAdjacentProfile(NULL)
{
}

//! Destructor
CGeomProfile::~CGeomProfile(void)
{
}

//! Returns this profile's coast ID
int CGeomProfile::nGetCoastID(void) const
{
   return m_nCoast;
}

//! Returns the profile's this-coast ID
int CGeomProfile::nGetProfileID(void) const
{
   return m_nProfileID;
}

//! Returns the coast point at which the profile starts
int CGeomProfile::nGetCoastPoint(void) const
{
   return m_nCoastPoint;
}

//! Returns a pointer to the location of the cell (grid CRS) on which the profile starts
CGeom2DIPoint* CGeomProfile::pPtiGetStartPoint(void)
{
   return &m_VCellInProfile.front();
}

//! Returns a pointer to the location of the cell (grid CRS) on which the profile ends
CGeom2DIPoint* CGeomProfile::pPtiGetEndPoint(void)
{
   return &m_VCellInProfile.back();
}

//! Sets a switch to indicate whether this is a start-of-coast profile
void CGeomProfile::SetStartOfCoast(bool const bFlag)
{
   m_bIsStartOfCoast = bFlag;
}

//! Returns the switch to indicate whether this is a start-of-coast profile
bool CGeomProfile::bIsStartOfCoast(void) const
{
   return m_bIsStartOfCoast;
}

//! Sets a switch to indicate whether this is an end-of-coast profile
void CGeomProfile::SetEndOfCoast(bool const bFlag)
{
   m_bIsEndOfCoast = bFlag;
}

//! Returns the switch to indicate whether this is an end-of-coast profile
bool CGeomProfile::bIsEndOfCoast(void) const
{
   return m_bIsEndOfCoast;
}

//! Returns true if this is a start-of-coast or an end-of-coast profile
bool CGeomProfile::bIsStartOrEndOfCoast(void) const
{
   if (m_bIsStartOfCoast || m_bIsEndOfCoast)
      return true;

   return false;
}

//! Flag this profile as having a CShore problem
void CGeomProfile::SetCShoreProblem(void)
{
   m_bCShoreProblem = true;
}

//! Returns true if this profile has a CShore problem
bool CGeomProfile::bHasCShoreProblem(void) const
{
   return m_bCShoreProblem;
}

//! Set the profile's status
void CGeomProfile::SetProfileStatus(int const nStatus)
{
   m_nProfileStatus = nStatus;
}

//! Gets the profile's status
int CGeomProfile::nGetProfileStatus(void) const
{
   return m_nProfileStatus;
}

//! Returns true if this is a problem-free profile
bool CGeomProfile::bProfileOK(void) const
{
   if (m_nProfileStatus == PROFILE_STATUS_OK)
      return true;

   return false;
}

//! Sets points (external CRS) in the profile. Note that only two points, the start and end point, are initially stored each profile
void CGeomProfile::SetPointsInProfile(vector<CGeom2DPoint> const* VNewPoints)
{
   CGeomMultiLine::m_VPtPoints = *VNewPoints;
}

//! Sets a single point (external CRS) in the profile, returns false if the point is out of range
bool CGeomProfile::bSetPointInProfile(int const nPoint, double const dNewX, double const dNewY)
{
   if (nPoint >= static_cast<int>(CGeomMultiLine::m_VPtPoints.size()))
      return false;

   CGeomMultiLine::m_VPtPoints[nPoint] = CGeom2DPoint(dNewX, dNewY);
   return true;
}

//! Appends a point (external CRS) to the profile
void CGeomProfile::AppendPointInProfile(double const dNewX, double const dNewY)
{
   CGeomMultiLine::m_VPtPoints.push_back(CGeom2DPoint(dNewX, dNewY));
}

//! Appends a point (external CRS) to the profile (overloaded version)
void CGeomProfile::AppendPointInProfile(CGeom2DPoint const* pPt)
{
   CGeomMultiLine::m_VPtPoints.push_back(*pPt);
}

//! Inserts an intersection (at a point specified in external CRS, with a line segment) into the profile
bool CGeomProfile::bInsertIntersection(double const dX, double const dY, int const nSeg)
{
   // Note no safety check to see if nSeg < nGetNumLineSegments()
   vector<CGeom2DPoint>::iterator it;
   it = CGeomMultiLine::m_VPtPoints.begin();

   // Do the insertion
   CGeomMultiLine::m_VPtPoints.insert(it + nSeg + 1, CGeom2DPoint(dX, dY));

   // Now insert a line segment in the associated multi-line, this will inherit the profile/line seg details from the preceding line segment
   CGeomMultiLine::InsertLineSegmentWithInheritance(nSeg);

   return true;
}

//! Truncates the profile's CGeomLine (external CRS points)
void CGeomProfile::TruncateProfile(int const nSize)
{
   CGeomMultiLine::m_VPtPoints.resize(nSize);
}

// void CGeomProfile::TruncateAndbSetPointInProfile(int const nPoint, double const dNewX, double const dNewY)
// {
// CGeomMultiLine::m_VPtPoints.resize(nPoint+1);
// CGeomMultiLine::m_VPtPoints[nPoint] = CGeom2DPoint(dNewX, dNewY);
// }

// void CGeomProfile::ShowProfile(void) const
// {
// for (int n = 0; n < CGeomMultiLine::m_VPtPoints.size(); n++)
// {
// cout << n << " [" << CGeomMultiLine::m_VPtPoints[n].dGetX() << "][" << CGeomMultiLine::m_VPtPoints[n].dGetY() << "]" << endl;
// }
// }

//! Returns the number of external CRS points in the profile (only two, initally; and always just two for grid-edge profiles). These points are stored even if the profile is later marked as invalid
int CGeomProfile::nGetProfileSize(void) const
{
   return static_cast<int>(CGeomMultiLine::m_VPtPoints.size());
}

//! Returns a single point (external CRS) from the profile
CGeom2DPoint* CGeomProfile::pPtGetPointInProfile(int const n)
{
   return &CGeomMultiLine::m_VPtPoints[n];
}

//! Returns a given external CRS point from the profile, and all points after this
vector<CGeom2DPoint> CGeomProfile::PtVGetThisPointAndAllAfter(int const nStart)
{
   return vector<CGeom2DPoint>(CGeomMultiLine::m_VPtPoints.begin() + nStart, CGeomMultiLine::m_VPtPoints.end());
}

//! Removes a line segment from the profile
// void CGeomProfile::RemoveLineSegment(int const nPoint)
// {
// m_VPtPoints.erase(CGeomMultiLine::m_VPtPoints.begin() + nPoint);
// CGeomMultiLine::RemoveLineSegment(nPoint);
// }

//! Queries the profile: is the given point (external CRS) a profile point?
bool CGeomProfile::bIsPointInProfile(double const dX, double const dY)
{
   CGeom2DPoint const Pt(dX, dY);
   auto it = find(CGeomMultiLine::m_VPtPoints.begin(), CGeomMultiLine::m_VPtPoints.end(), &Pt);

   if (it != CGeomMultiLine::m_VPtPoints.end())
      return true;
   else
      return false;
}

//! Queries the profile: is the given point (external CRS) a profile point? If so, then it also returns the number of the point in the profile
bool CGeomProfile::bIsPointInProfile(double const dX, double const dY, int& nPoint)
{
   CGeom2DPoint const Pt(dX, dY);
   auto it = find(CGeomMultiLine::m_VPtPoints.begin(), CGeomMultiLine::m_VPtPoints.end(), &Pt);

   if (it != CGeomMultiLine::m_VPtPoints.end())
   {
      // Found, so return true and set nPoint to be the index of the point which was found
      nPoint = static_cast<int>(it - CGeomMultiLine::m_VPtPoints.begin());
      return true;
   }

   else
      return false;
}

// int CGeomProfile::nFindInsertionLineSeg(double const dInsertX, double const dInsertY)
// {
// for (int n = 0; n < CGeomMultiLine::m_VPtPoints.back(); n++)
// {
// double
// dThisX = CGeomMultiLine::m_VPtPoints[n].dGetX(),
// dThisY = CGeomMultiLine::m_VPtPoints[n].dGetY(),
// dNextX = CGeomMultiLine::m_VPtPoints[n+1].dGetX(),
// dNextY = CGeomMultiLine::m_VPtPoints[n+1].dGetY();
//
// bool
// bBetweenX = false,
// bBetweenY = false;
//
// if (dNextX >= dThisX)
// {
//          // Ascending
// if ((dInsertX >= dThisX) && (dInsertX <= dNextX))
// bBetweenX = true;
// }
// else
// {
//          // Descending
// if ((dInsertX >= dNextX) && (dInsertX <= dThisX))
// bBetweenX = true;
// }
//
// if (dNextY >= dThisY)
// {
//          // Ascending
// if ((dInsertY >= dThisY) && (dInsertY <= dNextY))
// bBetweenY = true;
// }
// else
// {
//          // Descending
// if ((dInsertY >= dNextY) && (dInsertY <= dThisY))
// bBetweenY = true;
// }
//
// if (bBetweenX && bBetweenY)
// return n;
// }
//
// return -1;
// }

// void CGeomProfile::AppendPointShared(bool const bShared)
// {
// m_bVShared.push_back(bShared);
// }

// bool CGeomProfile::bPointShared(int const n) const
// {
//    // TODO 055 No check to see if n < size()
// return m_bVShared[n];
// }

//! Sets the up-coast adjacent profile
void CGeomProfile::SetUpCoastAdjacentProfile(CGeomProfile* pProfile)
{
   m_pUpCoastAdjacentProfile = pProfile;
}

// //! Returns the up-coast adjacent profile, returns NULL if there is no up-coast adjacent profile
// CGeomProfile* CGeomProfile::pGetUpCoastAdjacentProfile(void) const
// {
//    return m_pUpCoastAdjacentProfile;
// }

//! Sets the down-coast adjacent profile
void CGeomProfile::SetDownCoastAdjacentProfile(CGeomProfile* pProfile)
{
   m_pDownCoastAdjacentProfile = pProfile;
}

//! Returns the down-coast adjacent profile, returns NULL if there is no down-coast adjacent profile
CGeomProfile* CGeomProfile::pGetDownCoastAdjacentProfile(void) const
{
   return m_pDownCoastAdjacentProfile;
}

//! Appends a cell (grid CRS) to the profile
void CGeomProfile::AppendCellInProfile(CGeom2DIPoint const* pPti)
{
   m_VCellInProfile.push_back(*pPti);
}

//! Appends a cell (grid CRS) to the profile (overloaded version)
void CGeomProfile::AppendCellInProfile(int const nX, int const nY)
{
   m_VCellInProfile.push_back(CGeom2DIPoint(nX, nY));
}

//! Sets the profile's vector of cells (grid CRS)
void CGeomProfile::SetCellsInProfile(vector<CGeom2DIPoint> const* VNewPoints)
{
   m_VCellInProfile = *VNewPoints;
}

//! Returns all cells (grid CRS) in the profile
vector<CGeom2DIPoint>* CGeomProfile::pPtiVGetCellsInProfile(void)
{
   return &m_VCellInProfile;
}

//! Returns a single cell (grid CRS) in the profile
CGeom2DIPoint* CGeomProfile::pPtiGetCellInProfile(int const n)
{
   // TODO 055 No check to see if n < size()
   return &m_VCellInProfile[n];
}

//! Returns the last cell (grid CRS) in the profile
CGeom2DIPoint* CGeomProfile::pPtiGetLastCellInProfile(void)
{
   // In grid CRS
   return &m_VCellInProfile.back();
}

//! Returns the first cell (grid CRS) in the profile
CGeom2DIPoint* CGeomProfile::pPtiGetFirstCellInProfile(void)
{
   // In grid CRS
   return &m_VCellInProfile.front();
}

//! Returns the number of cells in the profile
int CGeomProfile::nGetNumCellsInProfile(void) const
{
   return static_cast<int>(m_VCellInProfile.size());
}

//! Returns the index of the cell on this profile which has a sea depth which is just less than a given depth. If every cell on the profile has a sea depth which is less than the given depth it returns INT_NODATA
int CGeomProfile::nGetCellGivenDepth(CGeomRasterGrid const* pGrid, double const dDepthIn)
{
   // int nIndex = INT_NODATA; // If not found, i.e. if every profile cell has sea depth less than dDepthIn
   int nIndex = static_cast<int>(m_VCellInProfile.size())-1;      // TODO bodge. Note needs to be size()-1, see e.g. beach_within_polygon.cpp line 73

   for (unsigned int n = 0; n < m_VCellInProfile.size(); n++)
   {
      int const nX = m_VCellInProfile[n].nGetX();
      int const nY = m_VCellInProfile[n].nGetY();

      double const dCellDepth = pGrid->m_Cell[nX][nY].dGetSeaDepth();

      if (dCellDepth >= dDepthIn)
      {
         nIndex = n;

         if (n > 0)
            nIndex = n - 1; // Grid CRS units

         break;
      }
   }

   return nIndex;
}

//! Sets the deep-water wave height for this profile
void CGeomProfile::SetProfileDeepWaterWaveHeight(double const dWaveHeight)
{
   m_dDeepWaterWaveHeight = dWaveHeight;
}

//! Returns the deep-water wave height for this profile
double CGeomProfile::dGetProfileDeepWaterWaveHeight(void) const
{
   return m_dDeepWaterWaveHeight;
}

//! Sets the deep-water wave orientation for this profile
void CGeomProfile::SetProfileDeepWaterWaveAngle(double const dWaveAngle)
{
   m_dDeepWaterWaveAngle = dWaveAngle;
}

//! Returns the deep-water wave orientation for this profile
double CGeomProfile::dGetProfileDeepWaterWaveAngle(void) const
{
   return m_dDeepWaterWaveAngle;
}

//! Sets the deep-water wave period for this profile
void CGeomProfile::SetProfileDeepWaterWavePeriod(double const dWavePeriod)
{
   m_dDeepWaterWavePeriod = dWavePeriod;
}

//! Returns the deep-water wave period for this profile
double CGeomProfile::dGetProfileDeepWaterWavePeriod(void) const
{
   return m_dDeepWaterWavePeriod;
}

//! Returns true if this is an intervention profile
bool CGeomProfile::bIsIntervention(void) const
{
   return m_bIntervention;
}

// //! Returns the index of a given cell in the vector of profile cells; returns INT_NODATA if not found
// int CGeomProfile::nGetIndexOfCellInProfile(int const nX, int const nY)
// {
//    for (unsigned int n = 0; n < m_VCellInProfile.size(); n++)
//    {
//       if ((m_VCellInProfile[n].nGetX() == nX) && (m_VCellInProfile[n].nGetY() == nY))
//          return n;
//    }
//
//    return INT_NODATA;
// }

