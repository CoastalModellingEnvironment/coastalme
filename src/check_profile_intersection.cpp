
/*!
   \file check_profile_intersection.cpp
   \brief Checks all coastline-normal profiles for intersection, and truncates or modifies those that intersect
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
#include <assert.h>

#include <cstdio>
#include <cmath>

#include <iostream>
using std::endl;
using std::ios;

#include <algorithm>
using std::find;

#include <utility>
using std::pair;

#include <random>
using std::normal_distribution;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2d_point.h"
#include "2di_point.h"

//===============================================================================================================================
//! For all coasts, checks all coastline-normal profiles for intersection, and modifies those that intersect
//===============================================================================================================================
void CSimulation::CheckAllProfilesForIntersection(void)
{
   LogStream << endl << m_ulIter << ": Checking for profile intersections" << endl;

   int const nCoastLines = static_cast<int>(m_VCoast.size());

   // Do once for every coastline object
   for (int nCoast = 0; nCoast < nCoastLines; nCoast++)
   {
      int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      bool bChanged = true;
      int nPass = -1;
      while (bChanged)
      {
         bChanged = false;
         nPass++;
         // LogStream << "***************************** nPass = " << nPass << endl;
         if (nPass >= MAX_ALONG_COAST_PASSES)
         {
            LogStream << m_ulIter << ":\t " << WARN << " reached maximum number (" << MAX_ALONG_COAST_PASSES << ") of along-coast passes for coast-normal intersection checks, abandoning check" << endl;
            break;
         }

         // Search in alternate directions, first down-coast (in the direction of increasing coast point numbers) then up-coast
         int nFirstProfileStartCoastPoint;
         for (int nFirstSearchDirection = DIRECTION_DOWNCOAST; nFirstSearchDirection <= DIRECTION_UPCOAST; nFirstSearchDirection++)
         {
            // Don't include coast-end profiles in the first profile search
            if (nFirstSearchDirection == DIRECTION_DOWNCOAST)
               nFirstProfileStartCoastPoint = 1;
            else
               nFirstProfileStartCoastPoint = nCoastSize - 2;

            // In the current search direction, look for "first profiles"
            for (int nFirstCoastPoint = nFirstProfileStartCoastPoint; ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? nFirstCoastPoint < nCoastSize : nFirstCoastPoint >= 0); ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? nFirstCoastPoint++ : nFirstCoastPoint--))
            {
               if (! m_VCoast[nCoast].bIsProfileAtCoastPoint(nFirstCoastPoint))
                  // No first profile at this coast point
                  continue;

               // There is a first profile at this coast point
               CGeomProfile* pFirstProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nFirstCoastPoint);
               // int const nFirstProfile = pFirstProfile->nGetProfileID();

               // Safety check: don't check this first profile if it is a start- or end-of-coast profile
               if (pFirstProfile->bIsStartOrEndOfCoast())
               {
                  // LogStream << m_ulIter << ":\t coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, nFirstCoastPoint = " << nFirstCoastPoint << " first profile " << pFirstProfile->nGetProfileID() << " is a start- or end-of-coast profile, do not check" << endl << endl;

                  continue;
               }

               // Don't check this first profile if it is has a problem
               if (! pFirstProfile->bProfileOK())
               {
                  // LogStream << m_ulIter << ":\t coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, nFirstCoastPoint = " << nFirstCoastPoint << " first profile " << pFirstProfile->nGetProfileID() << " is not OK, do not check" << endl;

                  continue;
               }

               // Does this profile have too many line segments?
               if (pFirstProfile->nGetNumLineSegments() > MAX_SEGMENTS)
               {
                  // Yes, too many line segments, so mark profile as invalid
                  pFirstProfile->SetProfileStatus(PROFILE_STATUS_TOO_MANY_SEGMENTS);

                  // LogStream << m_ulIter << ":\t coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, nFirstCoastPoint = " << nFirstCoastPoint << " first profile " << pFirstProfile->nGetProfileID() << " has more than " << MAX_SEGMENTS << " line segments, marking as invalid" << endl << endl;

                  continue;
               }

               // Now search in alternate directions, starting from one point either side of the first profile coast point, first down-coast (in the direction of increasing coast point numbers) then up-coast
               for (int nSecondSearchDirection = DIRECTION_DOWNCOAST; nSecondSearchDirection <= DIRECTION_UPCOAST; nSecondSearchDirection++)
               {
                  int nSecondProfileStartCoastPoint;

                  if (nSecondSearchDirection == DIRECTION_DOWNCOAST)
                     nSecondProfileStartCoastPoint = nFirstCoastPoint + 1;
                  else
                     nSecondProfileStartCoastPoint = nFirstCoastPoint - 1;

                  // Search for a "second profile"
                  for (int nSecondCoastPoint = nSecondProfileStartCoastPoint; (nSecondSearchDirection == DIRECTION_DOWNCOAST) ? nSecondCoastPoint < nCoastSize : nSecondCoastPoint >= 0; (nSecondSearchDirection == DIRECTION_DOWNCOAST) ? nSecondCoastPoint++ : nSecondCoastPoint--)
                  {
                     if (! m_VCoast[nCoast].bIsProfileAtCoastPoint(nSecondCoastPoint))
                        // No second profile at this coast point
                        continue;

                     // There is a second profile at this coast point, so get a pointer to the profile
                     CGeomProfile* pSecondProfile = m_VCoast[nCoast].pGetProfileAtCoastPoint(nSecondCoastPoint);
                     // int const nSecondProfile = pSecondProfile->nGetProfileID();

                     // LogStream << m_ulIter << ":\t  coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, first profile = " << pFirstProfile->nGetProfileID() << " (coast point " << nFirstCoastPoint << ") second profile = " << pSecondProfile->nGetProfileID() << " (coast point " << nSecondCoastPoint << ")" << endl;

                     // Don't check this second profile if it is has a problem
                     if (! pSecondProfile->bProfileOK())
                     {
                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, first profile = " << pFirstProfile->nGetProfileID() << " second profile " << pSecondProfile->nGetProfileID() << " is not OK, do not check" << endl;

                        continue;
                     }

                     // The first and second profiles have not previously intersected. So check for intersection now
                     int nProf1LineSeg = 0;
                     int nProf2LineSeg = 0;
                     double dIntersectX = 0;
                     double dIntersectY = 0;
                     double dAvgEndX = 0;
                     double dAvgEndY = 0;

                     if (! bCheckForIntersection(pFirstProfile, pSecondProfile, nProf1LineSeg, nProf2LineSeg, dIntersectX, dIntersectY, dAvgEndX, dAvgEndY))
                        // The first and second profiles do not intersect
                        continue;

                     // The first and second profiles do intersect. If the point of intersertion is present in both profiles then we must have we already dealt with this intersection
                     if ((pFirstProfile->bIsPointInProfile(dIntersectX, dIntersectY)) && (pSecondProfile->bIsPointInProfile(dIntersectX, dIntersectY)))
                        // Yes, present in both profiles
                        continue;

                     // The first and second profiles do intersect. Get the intersection point and the average endpoint in the grid CRS
                     // int nPoint = -1;
                     int nIntersectX = nRound(dExtCRSXToGridX(dIntersectX));
                     int nIntersectY = nRound(dExtCRSYToGridY(dIntersectY));
                     int nAvgEndX = nRound(dExtCRSXToGridX(dAvgEndX));
                     int nAvgEndY = nRound(dExtCRSYToGridY(dAvgEndY));

                     // Safety check: make sure that the point of intersection is within the valid grid
                     if (! bIsWithinValidGrid(nIntersectX, nIntersectY))
                     {
                        // LogStream << m_ulIter << ":\t  << coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, intersection {" << dIntersectX << ", " << dIntersectY << "} constrained to be ";

                        KeepWithinValidGrid(dIntersectX, dIntersectY, nIntersectX, nIntersectY);

                        // LogStream << "{" << dIntersectX << ", " << dIntersectY << "}" << endl;
                     }

                     // Safety check: make sure that the average endpoint is within the valid grid
                     if (! bIsWithinValidGrid(nAvgEndX, nAvgEndY))
                     {
                        // LogStream << m_ulIter << ":\t  << coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, average endpoint {" << dAvgEndX << ", " << dAvgEndY << "} constrained to be ";

                        KeepWithinValidGrid(dAvgEndX, dAvgEndY, nAvgEndX, nAvgEndY);

                        // LogStream << "{" << dAvgEndX << ", " << dAvgEndY << "}" << endl;
                     }

                     // Is the point of intersection also the end point of the first profile and the second profile? (Uncommon, but it happens occasionally)
                     int const nFirstPointSize = pFirstProfile->nGetProfileSize();
                     CGeom2DPoint const* pPtFirstEnd = pFirstProfile->pPtGetPointInProfile(nFirstPointSize-1);

                     int const nSecondPointSize = pSecondProfile->nGetProfileSize();
                     CGeom2DPoint const* pPtSecondEnd = pSecondProfile->pPtGetPointInProfile(nSecondPointSize-1);

                     if (bFPIsEqual(dIntersectX, pPtFirstEnd->dGetX(), TOLERANCE) && bFPIsEqual(dIntersectY, pPtFirstEnd->dGetY(), TOLERANCE) && bFPIsEqual(dIntersectX, pPtSecondEnd->dGetX(), TOLERANCE) && bFPIsEqual(dIntersectY, pPtSecondEnd->dGetY(), TOLERANCE))
                     {
                        // Yes, the point of intersection also the end point of the first profile and the second profile
                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, intersection point {" << dIntersectX << ", " << dIntersectY << "} already at seaward end of both profiles, no need to modify" << endl;

                        continue;
                     }

                     // The point of intersection is not also the endpoint of both profiles, so at least one profile must be modified
                     bChanged = true;

                     // Is the first profile an intervention profile?
                     if (pFirstProfile->bIsIntervention())
                     {
                        // Truncate the first profile, since it is an intervention profile
                        // LogStream << m_ulIter << ": profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " intersect, truncate first profile " << pFirstProfile->nGetProfileID() << " at {" << dIntersectX << ", " << dIntersectY << "} since it is an intervention profile" << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg);
#ifdef _DEBUG
                        // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
#endif
                        continue;
                     }

                     // Is the second profile an intervention profile?
                     if (pSecondProfile->bIsIntervention())
                     {
                        // Truncate the second profile, since it is an intervention profile
                        // LogStream << m_ulIter << ": profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " intersect, truncate second profile " << pSecondProfile->nGetProfileID() << " at {" << dIntersectX << ", " << dIntersectY << "} since it is an intervention profile" << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg);
#ifdef _DEBUG
                        // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
#endif
                        continue;
                     }

                     // Is the second profile a stsrt-of-cost or end-of-coast profile?
                     if (pSecondProfile->bIsStartOrEndOfCoast())
                     {
                        // It is, so truncate the first profile
                        // LogStream << m_ulIter << ": profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " intersect, and profile " << pSecondProfile->nGetProfileID() << " is a start-of-coast or end-of-coast prifile, so truncate first profile " << pFirstProfile->nGetProfileID() << " at {" << dIntersectX << ", " << dIntersectY << "} since it is an intervention profile" << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg);
#ifdef _DEBUG
                        // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
#endif
                        continue;

                     }



                     // The point of intersection is not already present in either profile, so get the number of line segments of each profile
                     int const nFirstProfileLineSegments = pFirstProfile->nGetNumLineSegments();
                     int const nSecondProfileLineSegments = pSecondProfile->nGetNumLineSegments();

                     // // DEBUG CODE ====================================================
                     // if (nFirstProfileLineSegments > 20)
                     // {
                     //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                     //    LogStream << "FIRST PROFILE ************************" << endl;
                     // }
                     //
                     // if (nSecondProfileLineSegments > 20)
                     // {
                     //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                     //    LogStream << "SECOND PROFILE ************************" << endl;
                     // }
                     // // DEBUG CODE ====================================================

                     // Next check whether the point of intersection is in the final line segment of both profiles
                     if ((nProf1LineSeg == (nFirstProfileLineSegments-1)) && (nProf2LineSeg == (nSecondProfileLineSegments-1)))
                     {
                        // Yes, the point of intersection is on the final line segment of both profiles
                        // LogStream << m_ulIter << ":\t  coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, end-segment intersection between profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " at [" << nIntersectX << "][" << nIntersectY << "] = {" << dIntersectX << ", " << dIntersectY << "} in line segment [" << nProf1LineSeg << "] of " << nFirstProfileLineSegments << " and line segment [" << nProf2LineSeg << "] of " << nSecondProfileLineSegments << " respectively" << endl;

                        // Merge the profiles seaward of the point of intersection
                        MergeProfilesAtFinalLineSegments(nCoast, pFirstProfile, pSecondProfile, nFirstProfileLineSegments, nSecondProfileLineSegments, dIntersectX, dIntersectY, dAvgEndX, dAvgEndY);

                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, merged first profile " << pFirstProfile->nGetProfileID() << " second profile = " << pSecondProfile->nGetProfileID() << " at {" << dIntersectX << ", " << dIntersectY << "}, averaged endpoint is {" << dAvgEndX << ", " << dAvgEndY << "}" << endl;
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                        continue;
                     }

                     // The profiles intersect, but the point of intersection is not in the final line segment of both profiles. One of the profiles will be truncated, the other profile will be retained
                     // LogStream << m_ulIter << ":\t  coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, intersection (not in both end segments) between profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " at [" << nIntersectX << "][" << nIntersectY << "] = {" << dIntersectX << ", " << dIntersectY << "} in line segment [" << nProf1LineSeg << "] of " << nFirstProfileLineSegments << " and line segment [" << nProf2LineSeg << "] of " << nSecondProfileLineSegments << ", respectively" << endl;

                     // Decide which profile to truncate, and which to retain
                     if (pFirstProfile->bIsIntervention())
                     {
                        // Truncate the first profile, since it is an intervention profile
                        // LogStream << m_ulIter << ":\t  coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, first profile = " << pFirstProfile->nGetProfileID() << " second profile = " << pSecondProfile->nGetProfileID() << ", first profile is an intervention profile, so truncate first profile (" << pFirstProfile->nGetProfileID() << ") at {" << dIntersectX << ", " << dIntersectY << "}" << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg);
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                        continue;
                     }

                     if (pSecondProfile->bIsIntervention())
                     {
                        // Truncate the second profile, since it is an intervention profile
                        if ((nFirstProfileLineSegments > 20) || (nSecondProfileLineSegments > 20))
                           // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, first profile = " << pFirstProfile->nGetProfileID() << " second profile = " << pSecondProfile->nGetProfileID() << ", second profile is an intervention profile, so truncate second profile (" << pSecondProfile->nGetProfileID() << ") at {" << dIntersectX << ", " << dIntersectY << "}" << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg);
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                        continue;
                     }

                     if (nFirstProfileLineSegments < nSecondProfileLineSegments)
                     {
                        // Truncate the first profile, since it has a smaller number of line segments
                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, first profile = " << pFirstProfile->nGetProfileID() << " has a smaller number of line segments, so truncate profile " << pFirstProfile->nGetProfileID() << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg);
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                        continue;
                     }

                     if (nFirstProfileLineSegments > nSecondProfileLineSegments)
                     {
                        // Truncate the second profile, since it has a smaller number of line segments
                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, second profile " << pSecondProfile->nGetProfileID() << " has a smaller number of line segments, so truncate profile " << pSecondProfile->nGetProfileID() << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg);
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                        continue;
                     }

                     // Both profiles have the same number of line segments, so choose randomly. Draw a sample from the unit normal distribution using random number generator 1
                     double const dRand = m_dGetFromUnitNormalDist(m_Rand[1]);
                     if (dRand >= 0)
                     {
                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, dRand = " << dRand << ", same number of line segments so randomly truncate first profile = " << pFirstProfile->nGetProfileID() << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pFirstProfile, pSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg);
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                     }
                     else
                     {
                        // LogStream << m_ulIter << ":\t   coast = " << nCoast << " nPass = " << nPass << " " << ((nFirstSearchDirection == DIRECTION_DOWNCOAST) ? "down-" : "up-") << "coast first profile search, " << (nSecondSearchDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast second profile search, dRand = " << dRand << ", same number of line segments so randomly truncate second profile = " << pSecondProfile->nGetProfileID() << endl;

                        TruncateOneProfileRetainOtherProfile(nCoast, pSecondProfile, pFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg);
#ifdef _DEBUG
                        // // DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        // // DEBUG CODE ====================================================
                        // if (nFirstProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "FIRST PROFILE ************************" << endl;
                        // }
                        //
                        // if (nSecondProfileLineSegments > 20)
                        // {
                        //    DEBUG_PrintProfileDetails(pFirstProfile, pSecondProfile);
                        //    LogStream << "SECOND PROFILE ************************" << endl;
                        // }
                        // // DEBUG CODE ====================================================
#endif
                     }
                  }
               }
            }
         }
      }
   }
}

//===============================================================================================================================
//! For all coasts, do further checks on all coastline-normal profiles
//===============================================================================================================================
int CSimulation::nFurtherCheckAndMarkAllProfiles(void)
{
   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      for (int n = 0; n < m_VCoast[nCoast].nGetNumProfiles(); n++)
      {
         CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
         int const nProfile = pProfile->nGetProfileID();

         if (pProfile->bProfileOK())
         {
            int const nSize = pProfile->nGetProfileSize();

            // Safety check
            if (nSize == 0)
            {
               pProfile->SetProfileStatus(PROFILE_STATUS_TOO_SHORT);

               LogStream << "Profile " << nProfile << " has zero length" << endl;
               continue;
            }

            CGeom2DPoint const* pPtEnd = pProfile->pPtGetPointInProfile(nSize - 1);
            CGeom2DIPoint const PtiEnd = PtiExtCRSToGridRound(pPtEnd);
            int nXEnd = PtiEnd.nGetX();
            int nYEnd = PtiEnd.nGetY();

            // Safety checks: the point may be outside the grid, so keep it within the grid
            nXEnd = tMin(nXEnd, m_nXGridSize - 1);
            nYEnd = tMin(nYEnd, m_nYGridSize - 1);
            nXEnd = tMax(nXEnd, 0);
            nYEnd = tMax(nYEnd, 0);

            // Is the water depth at the end point less than the depth of closure? We do this again because some profiles may have been shortened as a result of intersection. Do once for every coastline object
            if (m_pRasterGrid->m_Cell[nXEnd][nYEnd].dGetSeaDepth() < m_dDepthOfClosure)
            {
               pProfile->SetProfileStatus(PROFILE_STATUS_TOO_SHORT);

               // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               //    LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, is too short for depth of closure " << m_dDepthOfClosure << " at end point [" << nXEnd << "][" << nYEnd << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}, flagging as too short" << endl;
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Checks all line segments of a pair of coastline-normal profiles for intersection. If the lines intersect, returns true with the numbers of the line segments at which intersection occurs in nProfile1LineSegment and nProfile1LineSegment, the intersection point (external CRS) in dXIntersect and dYIntersect, and the 'average' seaward endpoint of the two intersecting profiles (external CRS) in dXAvgEnd and dYAvgEnd
//===============================================================================================================================
bool CSimulation::bCheckForIntersection(CGeomProfile* const pVProfile1, CGeomProfile* const pVProfile2, int& nProfile1LineSegment, int& nProfile2LineSegment, double& dXIntersect, double& dYIntersect, double& dXAvgEnd, double& dYAvgEnd)
{
   // For both profiles, look at all line segments
   int const nProfile1NumSegments = pVProfile1->nGetNumLineSegments();
   int const nProfile2NumSegments = pVProfile2->nGetNumLineSegments();

   for (int i = 0; i < nProfile1NumSegments; i++)
   {
      for (int j = 0; j < nProfile2NumSegments; j++)
      {
         // In external coordinates
         double const dX1 = pVProfile1->pPtVGetPoints()->at(i).dGetX();
         double const dY1 = pVProfile1->pPtVGetPoints()->at(i).dGetY();
         double const dX2 = pVProfile1->pPtVGetPoints()->at(i + 1).dGetX();
         double const dY2 = pVProfile1->pPtVGetPoints()->at(i + 1).dGetY();

         double const dX3 = pVProfile2->pPtVGetPoints()->at(j).dGetX();
         double const dY3 = pVProfile2->pPtVGetPoints()->at(j).dGetY();
         double const dX4 = pVProfile2->pPtVGetPoints()->at(j + 1).dGetX();
         double const dY4 = pVProfile2->pPtVGetPoints()->at(j + 1).dGetY();

         // Uses Cramer's Rule to solve the equations. Modified from code at http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect (in turn based on Andre LeMothe's "Tricks of the Windows Game Programming Gurus")
         double const dDiffX1 = dX2 - dX1;
         double const dDiffY1 = dY2 - dY1;
         double const dDiffX2 = dX4 - dX3;
         double const dDiffY2 = dY4 - dY3;

         double dS = -999;
         double dT = -999;
         double const dTmp = -dDiffX2 * dDiffY1 + dDiffX1 * dDiffY2;

         if (! bFPIsEqual(dTmp, 0.0, TOLERANCE))
         {
            dS = (-dDiffY1 * (dX1 - dX3) + dDiffX1 * (dY1 - dY3)) / dTmp;
            dT = (dDiffX2 * (dY1 - dY3) - dDiffY2 * (dX1 - dX3)) / dTmp;
         }

         if (dS >= 0 && dS <= 1 && dT >= 0 && dT <= 1)
         {
            // Collision detected, calculate intersection coordinates
            dXIntersect = dX1 + (dT * dDiffX1);
            dYIntersect = dY1 + (dT * dDiffY1);

            // And calc the average end-point coordinates
            dXAvgEnd = (dX2 + dX4) / 2;
            dYAvgEnd = (dY2 + dY4) / 2;

            // Get the line segments at which intersection occurred
            nProfile1LineSegment = i;
            nProfile2LineSegment = j;

            // LogStream << "\t" << "INTERSECTION dX2 = " << dX2 << " dX4 = " << dX4 << " dY2 = " << dY2 << " dY4 = " << dY4 << endl;
            return true;
         }
      }
   }

   // No intersection
   return false;
}

//===============================================================================================================================
//! For all coastlines, marks all coastline-normal profiles (apart from the two 'special' ones at the start and end of the coast) onto the raster grid, i.e. rasterizes multi-line vector objects onto the raster grid. Note that this doesn't work if the vector has already been interpolated to fit on the grid i.e. if distances between vector points are just one cell apart
//===============================================================================================================================
int CSimulation::nMarkProfilesOnGrid(void)
{
   int nValidProfiles = 0;

   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      // How many profiles on this coast?
      int const nProfiles = m_VCoast[nCoast].nGetNumProfiles();

      if (nProfiles == 0)
      {
         // This can happen if the coastline is very short, so just give a warning and carry on with the next coastline
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << WARN << m_ulIter << ": coast " << nCoast << " has no profiles" << endl;

         continue;
      }

      static bool bDownCoast = true;

      // Now do this for every profile, alternate between up-coast and down-coast directions
      for (int n = 0; n < nProfiles; n++)
      {
         CGeomProfile* pProfile;

         if (bDownCoast)
            pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
         else
            pProfile = m_VCoast[nCoast].pGetProfileWithUpCoastSeq(n);

         // Don't do this for the first and last profiles (i.e. the profiles at the start and end of the coast) since these are put onto the grid elsewhere
         if (pProfile->bIsStartOrEndOfCoast())
            continue;

         int const nProfile = pProfile->nGetProfileID();

         // If this profile has a problem, then forget about it
         if (! pProfile->bProfileOK())
         {
            // LogStream << m_ulIter << ": in nMarkProfilesOnGrid(), coast " << nCoast << " profile " << nProfile << " is not OK, not marked on grid" << endl;
            continue;
         }

         int const nPoints = pProfile->nGetProfileSize();

         if (nPoints < 2)
         {
            // Need at least two points in the profile, so this profile is invalid: mark it
            pProfile->SetProfileStatus(PROFILE_STATUS_TOO_SHORT);

            if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, has only " << nPoints << " points" << endl;

            continue;
         }

         // OK, go for it: set up temporary vectors to hold the x-y coords (in grid CRS) of the cells which we will mark
         vector<CGeom2DIPoint> VCellsToMark;
         vector<bool> bVShared;
         bool bTooShort = false;
         // bool const bTruncatedSameCoast = false;
         bool bHitCoast = false;
         bool bHitLand = false;
         bool bHitIntervention = false;
         // bool bHitAnotherProfile = false;

         CreateRasterizedProfile(nCoast, pProfile, &VCellsToMark, &bVShared, bTooShort, /*bTruncatedSameCoast,*/ bHitCoast, bHitLand, bHitIntervention/*, bHitAnotherProfile*/);

         if ((/*bTruncatedSameCoast &&*/ (! ACCEPT_TRUNCATED_PROFILES)) || bTooShort || bHitCoast || bHitLand || bHitIntervention /*|| bHitAnotherProfile*/ || VCellsToMark.size() == 0)
         {
            VCellsToMark.clear();
            bVShared.clear();
            continue;
         }

         // This profile is fine
         nValidProfiles++;

         for (unsigned int k = 0; k < VCellsToMark.size(); k++)
         {
            // Ignore duplicate points
            if ((k > 0) && (VCellsToMark[k] == m_VCoast[nCoast].pGetProfile(nProfile)->pPtiGetLastCellInProfile()))
               continue;

            // Mark each cell in the raster grid
            int const nXTmp = VCellsToMark[k].nGetX();
            int const nYTmp = VCellsToMark[k].nGetY();
            m_pRasterGrid->m_Cell[nXTmp][nYTmp].SetCoastAndProfileID(nCoast, nProfile);

            // Store the raster grid coordinates in the profile object
            m_VCoast[nCoast].pGetProfile(nProfile)->AppendCellInProfile(nXTmp, nYTmp);
         }

         // Get the deep water wave height and orientation values at the end of the profile
         double const dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWaveHeight();
         double const dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWaveAngle();
         double const dDeepWaterWavePeriod = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWavePeriod();

         // And store them for this profile
         m_VCoast[nCoast].pGetProfile(nProfile)->SetProfileDeepWaterWaveHeight(dDeepWaterWaveHeight);
         m_VCoast[nCoast].pGetProfile(nProfile)->SetProfileDeepWaterWaveAngle(dDeepWaterWaveAngle);
         m_VCoast[nCoast].pGetProfile(nProfile)->SetProfileDeepWaterWavePeriod(dDeepWaterWavePeriod);
      }

      bDownCoast = ! bDownCoast;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Given a pointer to a coastline-normal profile, returns an output vector of cells which are 'under' every line segment of the profile. If there is a problem with the profile (e.g. a rasterized cell is dry land or coast, or the profile has to be truncated) then we pass this back as an error code
//===============================================================================================================================
void CSimulation::CreateRasterizedProfile(int const nCoast, CGeomProfile* pProfile, vector<CGeom2DIPoint>* pVIPointsOut, vector<bool>* pbVShared, bool& bTooShort, /*bool const& bTruncatedSameCoast,*/ bool& bHitCoast, bool& bHitLand, bool& bHitIntervention/*, bool& bHitAnotherProfile*/)
{
   int const nProfile = pProfile->nGetProfileID();
   int nSeg = 0;
   int const nNumSegments = pProfile->nGetNumLineSegments();

   // LogStream << m_ulIter << ": in CreateRasterizedProfile() *pPtiStart for profile " << nProfile << " is [" << pPtiStart->nGetX() << "][" << pPtiStart->nGetY() << "]" << endl;
   int nXStartLast = INT_NODATA;
   int nYStartLast = INT_NODATA;
   int nXEndLast = INT_NODATA;
   int nYEndLast = INT_NODATA;

   // Do for every segment of this profile
   for (nSeg = 0; nSeg < nNumSegments; nSeg++)
   {
      // Do once for every line segment
      CGeom2DIPoint PtiSegStart;

      if (nSeg == 0)
      {
         // If this is the first segment, use the coastline start point to prevent external CRS to grid CRS rounding errors
         int const nCoastPoint = pProfile->nGetCoastPoint();
         PtiSegStart = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint);
      }
      else
      {
         CGeom2DPoint const* pPtSegStart = pProfile->pPtGetPointInProfile(nSeg);

         // Convert from the external CRS to grid CRS
         PtiSegStart = PtiExtCRSToGridRound(pPtSegStart);
      }

      CGeom2DPoint const* pPtSegEnd = pProfile->pPtGetPointInProfile(nSeg + 1); // This is OK

      // Convert from the external CRS to grid CRS
      CGeom2DIPoint const PtiSegEnd = PtiExtCRSToGridRound(pPtSegEnd);

      // Safety check
      if (PtiSegStart == PtiSegEnd)
         continue;

      int const nXStart = PtiSegStart.nGetX();
      int const nYStart = PtiSegStart.nGetY();
      int const nXEnd = PtiSegEnd.nGetX();
      int const nYEnd = PtiSegEnd.nGetY();

      bool bShared = false;

      if (pProfile->nGetNumCoincidentPairsForLineSegment(nSeg) > 1)
      {
         bShared = true;

         // If this is the second or more of several coincident line segments (i.e. it has the same start and end points as the previous line segment) then ignore it
         if ((nXStart == nXStartLast) && (nYStart == nYStartLast) && (nXEnd == nXEndLast) && (nYEnd == nYEndLast))
            continue;
      }

      // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
      double dXInc = nXEnd - nXStart;
      double dYInc = nYEnd - nYStart;
      double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

      dXInc /= dLength;
      dYInc /= dLength;

      double dX = nXStart;
      double dY = nYStart;

      // Process each interpolated point
      for (int m = 0; m <= nRound(dLength); m++)
      {
         int const nX = nRound(dX);
         int const nY = nRound(dY);

         // Safety check
         if (! bIsWithinValidGrid(nX, nY))
            continue;

         // Do some checking of this interpolated point, but only if this is not a grid-edge profile (these profiles are always valid)
         if (! pProfile->bIsStartOrEndOfCoast())
         {
            // Does this profile cross another profile (despite our best efforts so far to make sure that this doesn't happen)? Test by finding out whether this cell (or an adjacent cell: does not matter which) is already marked as 'under' a profile
            int nYTmp = nY+1;
            if (nY+1 >= m_nYGridSize)
               nYTmp = nY-1;

            // If this is the first line segment of the profile, then once we are clear of the coastline (when m > PROFILE_CHECK_DIST_FROM_COAST), check if this profile hits an intervention, or land, or another profile at this interpolated point. NOTE Get problems here since if the coastline vector has been heavily smoothed, this can result is 'false positives' profiles marked as invalid which are not actually invalid, because the profile hits land when m = 0 or m = 1. This results in some cells being flagged as profile cells which are actually inland
            if (m > PROFILE_CHECK_DIST_FROM_COAST)
            {
               // First check whether [nY][nY] hits an intervention
               if (m_pRasterGrid->m_Cell[nX][nY].nGetInterventionClass() != INT_NODATA)
               {
                  // We've hit an intervention, so set a switch and mark the profile
                  bHitIntervention = true;
                  pProfile->SetProfileStatus(PROFILE_STATUS_HIT_INTERVENTION);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " HIT INTERVENTION at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << ", SWL = " << m_dThisIterSWL << endl;

                  return;
               }

               // Next check whether [nY][nYTmp] hits an intervention
               if (m_pRasterGrid->m_Cell[nX][nYTmp].nGetInterventionClass() != INT_NODATA)
               {
                  // We've hit an intervention, so set a switch and mark the profile
                  bHitIntervention = true;
                  pProfile->SetProfileStatus(PROFILE_STATUS_HIT_INTERVENTION);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " HIT INTERVENTION at [" << nX << "][" << nYTmp << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nYTmp) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nYTmp].dGetAllSedTopElevIncTalus() << ", SWL = " << m_dThisIterSWL << endl;

                  return;
               }

               // Check whether [nY][nY] hits a coastline
               if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
               {
                  // We've hit a coastline so set a switch and mark the profile, then quit
                  bHitCoast = true;
                  pProfile->SetProfileStatus(PROFILE_STATUS_HIT_COAST);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " is invalid, HIT COAST " << m_pRasterGrid->m_Cell[nX][nY].nGetCoastline() << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

                  return;
               }

               // Check whether [nY][nYTmp] hits a coastline
               if (m_pRasterGrid->m_Cell[nX][nYTmp].bIsCoastline())
               {
                  // We've hit a coastline so set a switch and mark the profile, then quit
                  bHitCoast = true;
                  pProfile->SetProfileStatus(PROFILE_STATUS_HIT_COAST);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " is invalid, HIT COAST " << m_pRasterGrid->m_Cell[nX][nYTmp].nGetCoastline() << " at [" << nX << "][" << nYTmp << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nYTmp) << "}" << endl;

                  return;
               }

               // Check whether [nY][nY] hits land (such as an island)
               if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
               {
                  // We've hit dry land, so set a switch and mark the profile
                  bHitLand = true;
                  pProfile->SetProfileStatus(PROFILE_STATUS_HIT_LAND);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " HIT LAND at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() << ", SWL = " << m_dThisIterSWL << endl;

                  return;
               }

               // Check whether [nY][nYTmp] hits land (such as an island)
               if (! m_pRasterGrid->m_Cell[nX][nYTmp].bIsInContiguousSea())
               {
                  // We've hit dry land, so set a switch and mark the profile
                  bHitLand = true;
                  pProfile->SetProfileStatus(PROFILE_STATUS_HIT_LAND);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " HIT LAND at [" << nX << "][" << nYTmp << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nYTmp) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nYTmp].dGetAllSedTopElevIncTalus() << ", SWL = " << m_dThisIterSWL << endl;

                  return;
               }

               // Check to see if [nX][nY] hits another profile which is not a coincident normal to this normal
               if (m_pRasterGrid->m_Cell[nX][nY].bIsProfile())
               {
                  // We've hit a raster cell which is already marked as 'under' a normal profile. Get the number of the profile which marked this cell, and the coast to which this profile belongs
                  int const nHitProfile = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();
                  int const nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nY].nGetProfileCoastID();

                  // Do both profiles belong to the same coast?
                  if (nCoast == nHitProfileCoast)
                  {
                     // Both profiles belong to the same coast. Is this the number of a coincident profile of this profile?
                     if (! pProfile->bFindFirstFromCoincidentPairsInLastLineSegment(nHitProfile))
                     {
                        // It isn't a coincident profile, so we have just hit an unrelated profile. Mark this profile as invalid and move on
                        // pProfile->SetProfileStatus(PROFILE_STATUS_HIT_PROFILE);
                        // bHitAnotherProfile = true;
                        //
                        // if (m_nLogFileDetail >= LOG_FILE_ALL)
                        //    LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " is invalid, hit another profile (" << nHitProfile << ") belonging to coast " << nHitProfileCoast << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
                        //
                        // return;
                     }
                  }
               }

               // Check to see if [nX][nTmp] hits another profile which is not a coincident normal to this normal
               if (m_pRasterGrid->m_Cell[nX][nYTmp].bIsProfile())
               {
                  // We've hit a raster cell which is already marked as 'under' a normal profile. Get the number of the profile which marked this cell, and the coast to which this profile belongs
                  int const nHitProfile = m_pRasterGrid->m_Cell[nX][nYTmp].nGetProfileID();
                  int const nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nYTmp].nGetProfileCoastID();

                  // Do both profiles belong to the same coast?
                  if (nCoast == nHitProfileCoast)
                  {
                     // Both profiles belong to the same coast. Is this the number of a coincident profile of this profile?
                     if (! pProfile->bFindFirstFromCoincidentPairsInLastLineSegment(nHitProfile))
                     {
                        // It isn't a coincident profile, so we have just hit an unrelated profile. Mark this profile as invalid and move on
                        // pProfile->SetProfileStatus(PROFILE_STATUS_HIT_PROFILE);
                        // bHitAnotherProfile = true;
                        //
                        // if (m_nLogFileDetail >= LOG_FILE_ALL)
                        //    LogStream << m_ulIter << ":\t coast " << nCoast << " profile " << nProfile << " is invalid, hit another profile (" << nHitProfile << ") belonging to coast " << nHitProfileCoast << " at [" << nX << "][" << nYTmp << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nYTmp) << "}" << endl;
                        //
                        // return;
                     }
                  }
               }
            }
         }

         // All OK, so append this point to the output vector
         pVIPointsOut->push_back(CGeom2DIPoint(nX, nY)); // Is in raster grid coordinates
         pbVShared->push_back(bShared);

         // And increment for next time
         dX += dXInc;
         dY += dYInc;
      }

      nXStartLast = nXStart;
      nYStartLast = nYStart;
      nXEndLast = nXEnd;
      nYEndLast = nYEnd;

      // if (bTruncatedSameCoast)
      //    break;
   }

   // if (bTruncatedSameCoast)
   // {
   //    if (nSeg < (nNumSegments - 1))
   //       // We are truncating the profile, so remove any line segments after this one
   //       pProfile->TruncateLineSegments(nSeg);
   //
   //    // Shorten the vector input. Ignore CPPCheck errors here, since we know that pVIPointsOut is not empty
   //    int const nLastX = pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX();
   //    int const nLastY = pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY();
   //
   //    pProfile->pPtGetPointInProfile(nSeg + 1)->SetX(dGridCentroidXToExtCRSX(nLastX));
   //    pProfile->pPtGetPointInProfile(nSeg + 1)->SetY(dGridCentroidYToExtCRSY(nLastY));
   // }

   if (pVIPointsOut->size() < 3)
   {
      // Coastline-normal profiles cannot be very short (e.g. with less than 3 cells), since we cannot calculate along-profile slope properly for such short profiles
      bTooShort = true;
      pProfile->SetProfileStatus(PROFILE_STATUS_TOO_SHORT);

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         // Ignore CPPCheck errors here, since we know that pVIPointsOut is not empty
         LogStream << m_ulIter << ": profile " << nProfile << " is invalid, is too short, only " << pVIPointsOut->size() << " points, HitLand?" << bHitLand << ". From [" << pVIPointsOut->at(0).nGetX() << "][" << pVIPointsOut->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pVIPointsOut->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pVIPointsOut->at(0).nGetY()) << "} to [" << pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX() << "][" << pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX()) << ", " << dGridCentroidYToExtCRSY(pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY()) << "}" << endl;
      }
   }
}

//===============================================================================================================================
//! Merges two profiles which intersect at their final (most seaward) line segments, seaward of their point of intersection (external CRS)
//===============================================================================================================================
void CSimulation::MergeProfilesAtFinalLineSegments(int const nCoast, CGeomProfile* pFirstProfile, CGeomProfile* pSecondProfile, int const nFirstProfileLineSegments, int const nSecondProfileLineSegments, double const dIntersectX, double const dIntersectY, double const dAvgEndX, double const dAvgEndY)
{
   // The point of intersection is on the final (most seaward) line segment of both profiles. Put together a vector of coincident profile numbers (with no duplicates) for both profiles
   int nCombinedLastSeg = 0;
   vector<pair<int, int>> prVCombinedProfilesCoincidentProfilesLastSeg;

   for (unsigned int n = 0; n < pFirstProfile->pprVGetCoincidentPairsForLineSegment(nFirstProfileLineSegments - 1)->size(); n++)
   {
      pair<int, int> prTmp;
      prTmp.first = pFirstProfile->pprVGetCoincidentPairsForLineSegment(nFirstProfileLineSegments - 1)->at(n).first;
      prTmp.second = pFirstProfile->pprVGetCoincidentPairsForLineSegment(nFirstProfileLineSegments - 1)->at(n).second;

      bool bFound = false;

      for (unsigned int m = 0; m < prVCombinedProfilesCoincidentProfilesLastSeg.size(); m++)
      {
         if (prVCombinedProfilesCoincidentProfilesLastSeg[m].first == prTmp.first)
         {
            bFound = true;
            break;
         }
      }

      if (! bFound)
      {
         prVCombinedProfilesCoincidentProfilesLastSeg.push_back(prTmp);
         nCombinedLastSeg++;
      }
   }

   for (unsigned int n = 0; n < pSecondProfile->pprVGetCoincidentPairsForLineSegment(nSecondProfileLineSegments - 1)->size(); n++)
   {
      pair<int, int> prTmp;
      prTmp.first = pSecondProfile->pprVGetCoincidentPairsForLineSegment(nSecondProfileLineSegments - 1)->at(n).first;
      prTmp.second = pSecondProfile->pprVGetCoincidentPairsForLineSegment(nSecondProfileLineSegments - 1)->at(n).second;

      bool bFound = false;

      for (unsigned int m = 0; m < prVCombinedProfilesCoincidentProfilesLastSeg.size(); m++)
      {
         if (prVCombinedProfilesCoincidentProfilesLastSeg[m].first == prTmp.first)
         {
            bFound = true;
            break;
         }
      }

      if (! bFound)
      {
         prVCombinedProfilesCoincidentProfilesLastSeg.push_back(prTmp);
         nCombinedLastSeg++;
      }
   }

   // Increment the number of each line segment
   for (int m = 0; m < nCombinedLastSeg; m++)
      prVCombinedProfilesCoincidentProfilesLastSeg[m].second++;

   vector<pair<int, int>> prVFirstProfileCoincidentProfilesLastSeg = *pFirstProfile->pprVGetCoincidentPairsForLineSegment(nFirstProfileLineSegments - 1);
   vector<pair<int, int>> prVSecondProfileCoincidentProfilesLastSeg = *pSecondProfile->pprVGetCoincidentPairsForLineSegment(nSecondProfileLineSegments - 1);
   int const nNumFirstProfileCoincidentProfilesLastSeg = static_cast<int>(prVFirstProfileCoincidentProfilesLastSeg.size());
   int const nNumSecondProfileCoincidentProfilesLastSeg = static_cast<int>(prVSecondProfileCoincidentProfilesLastSeg.size());

   // LogStream << m_ulIter << ":\t   coast " << nCoast << " profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " end-segment intersection at line segment " << nFirstProfileLineSegments-1 << " of " << nFirstProfileLineSegments << ", and line segment " << nSecondProfileLineSegments-1 << " of " << nSecondProfileLineSegments << ", respectively. Both truncated at {" << dIntersectX << ", " << dIntersectY << "}, combined profiles " << pFirstProfile->nGetProfileID() << " and " << pSecondProfile->nGetProfileID() << " extended to {" << dAvgEndX << ", " << dAvgEndY << "}" << endl;

   // Truncate the first profile, and all co-incident profiles, at the point of intersection
   for (int n = 0; n < nNumFirstProfileCoincidentProfilesLastSeg; n++)
   {
      int const nThisProfile = prVFirstProfileCoincidentProfilesLastSeg[n].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);
      int const nProfileLength = pThisProfile->nGetProfileSize();

      // This is the final line segment of the first 'main' profile. We are assuming that it is also the final line segment of all co-incident profiles. This is fine, altho' each profile may well have a different number of line segments landwards i.e. the number of the line segment may be different for each co-incident profile
      if (! pThisProfile->bSetPointInProfile(nProfileLength - 1, dIntersectX, dIntersectY))
         break;
   }

   // Truncate the second profile, and all co-incident profiles, at the point of intersection
   for (int n = 0; n < nNumSecondProfileCoincidentProfilesLastSeg; n++)
   {
      int const nThisProfile = prVSecondProfileCoincidentProfilesLastSeg[n].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);
      int const nProfileLength = pThisProfile->nGetProfileSize();

      // This is the final line segment of the second 'main' profile. We are assuming that it is also the final line segment of all co-incident profiles. This is fine, altho' each profile may well have a different number of line segments landwards i.e. the number of the line segment may be different for each co-incident profile
      if (! pThisProfile->bSetPointInProfile(nProfileLength - 1, dIntersectX, dIntersectY))
         break;
   }

   // Append a new straight line segment to the existing line segment(s) of the first profile, and to all co-incident profiles
   for (int nThisLineSeg = 0; nThisLineSeg < nNumFirstProfileCoincidentProfilesLastSeg; nThisLineSeg++)
   {
      int const nThisProfile = prVFirstProfileCoincidentProfilesLastSeg[nThisLineSeg].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Append a new point to this first profile
      pThisProfile->AppendPointInProfile(dAvgEndX, dAvgEndY);

      // Append an empty line segment to the first profile's CGeomMultiLine then add all profile numbers from the previous line segment to this new line segment
      pThisProfile->AppendLineSegment();
      for (int m = 0; m < nCombinedLastSeg; m++)
         pThisProfile->AppendPairToFinalLineSegment(prVCombinedProfilesCoincidentProfilesLastSeg[m]);

      assert(pThisProfile->nGetProfileSize() == (1 + pThisProfile->nGetNumLineSegments()));
   }

   // Append a new straight line segment to the existing line segment(s) of the second profile, and to all co-incident profiles
   for (int nThisLineSeg = 0; nThisLineSeg < nNumSecondProfileCoincidentProfilesLastSeg; nThisLineSeg++)
   {
      int const nThisProfile = prVSecondProfileCoincidentProfilesLastSeg[nThisLineSeg].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Append a new point to this second profile
      pThisProfile->AppendPointInProfile(dAvgEndX, dAvgEndY);

      // Append an empty line segment to the second profile's CGeomMultiLine then add all profile numbers from the previous line segment to this new line segment
      pThisProfile->AppendLineSegment();
      for (int m = 0; m < nCombinedLastSeg; m++)
         pThisProfile->AppendPairToFinalLineSegment(prVCombinedProfilesCoincidentProfilesLastSeg[m]);
   }
}

//===============================================================================================================================
//! Truncates one intersecting profile at the point of intersection (external CRS), and retains the other profile
//===============================================================================================================================
void CSimulation::TruncateOneProfileRetainOtherProfile(int const nCoast, CGeomProfile* pProfileToTruncate, CGeomProfile* pProfileToRetain, double dIntersectX, double dIntersectY, int nProfileToTruncateIntersectLineSeg, int nProfileToRetainIntersectLineSeg)
{
   CGeom2DPoint const PtIntersect(dIntersectX, dIntersectY);
   bool const bAlreadyPresent = pProfileToRetain->bIsPresent(&PtIntersect);

   // Insert the intersection point into the main retain-profile if it is not already in the profile, and do the same for all co-incident profiles of the main retain-profile. Also add details of the to-truncate profile (and all its coincident profiles) to every line segment of the main to-retain profile which is seaward of the point of intersection
   int const nRet = nInsertPointIntoProfilesIfNeededThenUpdate(nCoast, pProfileToRetain, dIntersectX, dIntersectY, nProfileToRetainIntersectLineSeg, pProfileToTruncate, nProfileToTruncateIntersectLineSeg, bAlreadyPresent);
   if (nRet != RTN_OK)
   {
      LogStream << m_ulIter << ": error in nInsertPointIntoProfilesIfNeededThenUpdate()" << endl;
      return;
   }

   // Get all profile points of the main retain-profile seawards from the intersection point, and do the same for the corresponding line segments (including coincident profiles). This also includes details of the main to-truncate profile (and all its coincident profiles)
   vector<CGeom2DPoint> PtVProfileLastPart;
   vector<vector<pair<int, int>>> prVLineSegLastPart;

   if (bAlreadyPresent)
   {
      PtVProfileLastPart = pProfileToRetain->PtVGetThisPointAndAllAfter(nProfileToRetainIntersectLineSeg);
      prVLineSegLastPart = pProfileToRetain->prVVGetAllLineSegmentsAfter(nProfileToRetainIntersectLineSeg);
   }
   else
   {
      PtVProfileLastPart = pProfileToRetain->PtVGetThisPointAndAllAfter(nProfileToRetainIntersectLineSeg + 1);
      prVLineSegLastPart = pProfileToRetain->prVVGetAllLineSegmentsAfter(nProfileToRetainIntersectLineSeg + 1);
   }

   // Truncate the truncate-profile at the point of intersection, and do the same for all its co-incident profiles. Then append the profile points of the main to-retain profile seaward from the intersection point, and do the same for the corresponding line segments (including coincident profiles)
   TruncateProfileAndAppendNew(nCoast, pProfileToTruncate, nProfileToTruncateIntersectLineSeg, &PtVProfileLastPart, &prVLineSegLastPart);
}

//===============================================================================================================================
//! Inserts an intersection point into the profile that is to be retained, if that point is not already present in the profile, then does the same for all co-incident profiles. Finally adds the numbers of the to-truncate profile (and all its coincident profiles) to the seaward line segments of the to-retain profile and all its coincident profiles
//===============================================================================================================================
int CSimulation::nInsertPointIntoProfilesIfNeededThenUpdate(int const nCoast, CGeomProfile* pProfileToRetain, double const dIntersectX, double const dIntersectY, int const nProfileToRetainIntersectLineSeg, CGeomProfile* pProfileToTruncate, int const nProfileToTruncateIntersectLineSeg, bool const bAlreadyPresent)
{
   int const nProfileToRetain = pProfileToRetain->nGetProfileID();

   // Get the index numbers of all coincident profiles for the 'main' to-retain profile for the line segment in which intersection occurs
   vector<pair<int, int>> prVCoincidentProfiles = *pProfileToRetain->pprVGetCoincidentPairsForLineSegment(nProfileToRetainIntersectLineSeg);
   int const nNumCoincident = static_cast<int>(prVCoincidentProfiles.size());
   vector<int> VnLineSegAfterIntersect(nNumCoincident, INT_NODATA); // The line segment after the point of intersection, for each co-incident profile

   // Do this for the main profile and all profiles which are co-incident for this line segment
   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      int const nThisProfile = prVCoincidentProfiles[nn].first;  // The number of this profile
      int const nThisLineSeg = prVCoincidentProfiles[nn].second; // The line segment of this profile

      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Is the intersection point already present in the to-retain profile?
      if (! bAlreadyPresent)
      {
         if (nThisProfile != nProfileToRetain)
            continue;

         // Safety check
         if (nThisLineSeg >= pThisProfile->nGetNumLineSegments())
         {
            LogStream << WARN << m_ulIter << ": *** nThisLineSeg = " << nThisLineSeg << " pThisProfile->nGetNumLineSegments() = " << pThisProfile->nGetNumLineSegments() << endl;
            continue;
         }

         // It is not already present, so insert it and also update the associated multi-line
         if (! pThisProfile->bInsertIntersection(dIntersectX, dIntersectY, nThisLineSeg))
         {
            // Error
            LogStream << WARN << m_ulIter << ": cannot insert intersection point {" << dIntersectX << ", " << dIntersectY << "} nThisLineSeg = " << nThisLineSeg << " into profile " << pThisProfile->nGetCoastID() << " (main profile) which has length " << pThisProfile->nGetProfileSize() << ", abandoning" << endl;
            LogStream << "\t ";
            LogStream << "Profile points for profile " << pThisProfile->nGetProfileID() << endl;
            for (int mmm = 0; mmm < pThisProfile->nGetProfileSize(); mmm++)
            {
               CGeom2DPoint const* pPtThis = pThisProfile->pPtGetPointInProfile(mmm);
               LogStream << "{" << pPtThis->dGetX() << ", " << pPtThis->dGetY() << "}\t";
            }
            LogStream << endl;

            LogStream << "Line segments for profile " << pThisProfile->nGetProfileID() << endl;
            int const nNumLineSeg = pThisProfile->nGetNumLineSegments();
            for (int nSeg = 0; nSeg < nNumLineSeg; nSeg++)
            {
               LogStream << "Line segment " << nSeg << endl;
               vector<pair<int, int>> const VprTmp = *pThisProfile->pprVGetCoincidentPairsForLineSegment(nSeg);
               int const nNumPairedCoinc = static_cast<int>(VprTmp.size());

               for (int nCoinc = 0; nCoinc < nNumPairedCoinc; nCoinc++)
               {
                  pair<int, int> prTmp;
                  prTmp.first = pThisProfile->pprVGetCoincidentPairsForLineSegment(nSeg)->at(nCoinc).first;
                  prTmp.second = pThisProfile->pprVGetCoincidentPairsForLineSegment(nSeg)->at(nCoinc).second;

                  LogStream << "{" << prTmp.first << ", " << prTmp.second << "}\t" << endl;
               }
               LogStream << endl;
            }
            LogStream << endl;

            // return RTN_ERR_CANNOT_INSERT_POINT;
            continue;
         }
      }

      // Get the line segment after intersection
      VnLineSegAfterIntersect[nn] = nThisLineSeg + 1;
   }

   // Get the coincident profiles for the to-truncate profile, at the line segment where intersection occurs
   vector<pair<int, int>> prVToTruncateCoincidentProfiles = *pProfileToTruncate->pprVGetCoincidentPairsForLineSegment(nProfileToTruncateIntersectLineSeg);
   int const nNumToTruncateCoincident = static_cast<int>(prVToTruncateCoincidentProfiles.size());

   // Now add the number of the to-truncate profile, and all its coincident profiles, to all line segments which are seaward of the point of intersection. Do this for the main profile and all profiles which are co-incident for this line segment
   for (int nCoinc = 0; nCoinc < nNumCoincident; nCoinc++)
   {
      int const nThisProfile = prVCoincidentProfiles[nCoinc].first; // The number of this profile

      if (nThisProfile != nProfileToRetain)
         continue;

      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Get the number of line segments for this to-retain profile (will have just increased, if we just inserted a point)
      int const nNumLineSegs = pThisProfile->nGetNumLineSegments();

      // Do for all line segments seaward of the point of intersection
      for (int nLineSeg = VnLineSegAfterIntersect[nCoinc], nIncr = 0; nLineSeg < nNumLineSegs; nLineSeg++, nIncr++)
      {
         // Safety check TODO ugly
         if (nLineSeg == INT_NODATA)
         {
            LogStream << "*** NO_DATA" << endl;
            // break;
            return RTN_ERR_CHECK_INTERSECTION;
         }

         // Add the number of the to-truncate profile, and all its coincident profiles, to this line segment
         for (int m = 0; m < nNumToTruncateCoincident; m++)
         {
            int const nProfileToAdd = prVToTruncateCoincidentProfiles[m].first;
            int const nProfileToAddLineSeg = prVToTruncateCoincidentProfiles[m].second;

            pThisProfile->AddCoincidentPairToExistingLineSegmentIfNotAlready(nLineSeg, nProfileToAdd, nProfileToAddLineSeg + nIncr);
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Truncate a profile at the point of intersection, and do the same for all its co-incident profiles
//===============================================================================================================================
void CSimulation::TruncateProfileAndAppendNew(int const nCoast, CGeomProfile* pProfileToTruncate, int const nMainProfileIntersectLineSeg, vector<CGeom2DPoint> const* pPtVProfileLastPart, vector<vector<pair<int, int>>> const* pprVLineSegLastPart)
{
   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurs
   vector<pair<int, int>> prVCoincidentProfiles = *pProfileToTruncate->pprVGetCoincidentPairsForLineSegment(nMainProfileIntersectLineSeg);
   int const nNumCoincident = static_cast<int>(prVCoincidentProfiles.size());

   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      // Do this for the main to-truncate profile, and do the same for all its co-incident profiles
      int const nThisProfile = prVCoincidentProfiles[nn].first;
      int const nThisProfileLineSeg = prVCoincidentProfiles[nn].second;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Truncate the profile
      // LogStream << "\tTruncating " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} at line segment " << nThisProfileLineSeg+1 << endl;
      pThisProfile->TruncateProfile(nThisProfileLineSeg + 1);

      // Reduce the number of line segments for this profile
      pThisProfile->TruncateLineSegments(nThisProfileLineSeg + 1);

      // Append the profile points from the last part of the retain-profile
      for (unsigned int mm = 0; mm < pPtVProfileLastPart->size(); mm++)
      {
         CGeom2DPoint const Pt = pPtVProfileLastPart->at(mm);
         pThisProfile->AppendPointInProfile(&Pt);
      }

      // Append to this profile all line segments, and their co-incident profile numbers, from the last part of the retain-profile
      for (unsigned int mm = 0; mm < pprVLineSegLastPart->size(); mm++)
      {
         vector<pair<int, int>> prVTmp = pprVLineSegLastPart->at(mm);

         pThisProfile->AppendLineSegment(&prVTmp);
      }

      // Fix the line segment numbers for this profile
      vector<int> nVProf;
      vector<int> nVProfsLineSeg;

      for (int nSeg = 0; nSeg < pThisProfile->nGetNumLineSegments(); nSeg++)
      {
         for (int nCoinc = 0; nCoinc < pThisProfile->nGetNumCoincidentPairsForLineSegment(nSeg); nCoinc++)
         {
            int const nProf = pThisProfile->nGetPairFirst(nSeg, nCoinc);
            int const nProfsLineSeg = pThisProfile->nGetPairSecond(nSeg, nCoinc);

            auto it = find(nVProf.begin(), nVProf.end(), nProf);

            if (it == nVProf.end())
            {
               // Not found
               nVProf.push_back(nProf);
               nVProfsLineSeg.push_back(nProfsLineSeg);
            }
            else
            {
               // Found
               int const nPos = static_cast<int>(it - nVProf.begin());
               int nNewProfsLineSeg = nVProfsLineSeg[nPos];
               nNewProfsLineSeg++;

               nVProfsLineSeg[nPos] = nNewProfsLineSeg;
               pThisProfile->SetPairSecond(nSeg, nCoinc, nNewProfsLineSeg);
            }
         }
      }
      assert(pThisProfile->nGetProfileSize() == (1 + pThisProfile->nGetNumLineSegments()));
   }
}

// #ifdef _DEBUG
// //===============================================================================================================================
// //! DEBUG ONLY: print profile details to logfile
// //===============================================================================================================================
// void CSimulation::DEBUG_PrintProfileDetails(CGeomProfile* pFirstProfile, CGeomProfile* pSecondProfile)
// {
//    // if ((nFirstProfileLineSegments > 20) || (nSecondProfileLineSegments > 20))
//    {
//       m_nExtra++;
//       LogStream << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% m_ulIter = " << m_ulIter << " m_nExtra = " << m_nExtra << endl;
//       LogStream << "first profile = " << pFirstProfile->nGetProfileID() << " second profile = " << pSecondProfile->nGetProfileID() << endl;
//
//       int const nFirstProfilePoints = pFirstProfile->nGetProfileSize();
//       int const nSecondProfilePoints = pSecondProfile->nGetProfileSize();
//       int const nFirstProfileLineSeg = pFirstProfile->nGetNumLineSegments();
//       int const nSecondProfileLineSeg = pSecondProfile->nGetNumLineSegments();
//
//       LogStream << "\tFirst profile = " << pFirstProfile->nGetProfileID() << " now has " << nFirstProfilePoints << " points" << endl;
//       LogStream << "\tPoints for profile " << pFirstProfile->nGetProfileID() << " are ";
//       for (int m = 0; m < nFirstProfilePoints; m++)
//          LogStream << "{" << pFirstProfile->pPtGetPointInProfile(m)->dGetX() << ", " << pFirstProfile->pPtGetPointInProfile(m)->dGetY() << "} ";
//       LogStream << endl;
//
//       LogStream << "\tFirst profile = " << pFirstProfile->nGetProfileID() << " now has " << nFirstProfileLineSeg << " line segments" << endl;
//       for (int m = 0; m < nFirstProfileLineSeg; m++)
//       {
//          vector<pair<int, int> > prVCoincidentProfiles = *pFirstProfile->pprVGetCoincidentPairsForLineSegment(m);
//          LogStream << "\tCo-incident profiles and line segments for line segment " << m << " of profile " << pFirstProfile->nGetProfileID() << " are ";
//          for (int nn = 0; nn < static_cast<int>(prVCoincidentProfiles.size()); nn++)
//             LogStream << "{" << prVCoincidentProfiles[nn].first << ", " << prVCoincidentProfiles[nn].second << "} ";
//          LogStream << " " << endl;
//       }
//
//       LogStream << "\tSecond profile = " << pSecondProfile->nGetProfileID() << " now has " << nSecondProfilePoints << " points" << endl;
//       LogStream << "\tPoints for profile " << pSecondProfile->nGetProfileID() << " are ";
//       for (int m = 0; m < nSecondProfilePoints; m++)
//          LogStream << "{" << pSecondProfile->pPtGetPointInProfile(m)->dGetX() << ", " << pSecondProfile->pPtGetPointInProfile(m)->dGetY() << "} ";
//       LogStream << endl;
//
//       LogStream << "\tSecond profile = " << pSecondProfile->nGetProfileID() << " now has " << nSecondProfileLineSeg << " line segments" << endl;
//       for (int m = 0; m < nSecondProfileLineSeg; m++)
//       {
//          vector<pair<int, int> > prVCoincidentProfiles = *pSecondProfile->pprVGetCoincidentPairsForLineSegment(m);
//          LogStream << "\tCo-incident profiles and line segments for line segment " << m << " of profile " << pSecondProfile->nGetProfileID() << " are ";
//          for (int nn = 0; nn < static_cast<int>(prVCoincidentProfiles.size()); nn++)
//             LogStream << "{" << prVCoincidentProfiles[nn].first << ", " << prVCoincidentProfiles[nn].second << "} ";
//          LogStream << " " << endl;
//       }
//
//       LogStream << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
//
//       string const strExtra = "_" + to_string(m_nExtra);
//       bWriteVectorGISFile(VECTOR_PLOT_NORMALS, &VECTOR_PLOT_NORMALS_TITLE, strExtra);
//       bWriteVectorGISFile(VECTOR_PLOT_INVALID_NORMALS, &VECTOR_PLOT_INVALID_NORMALS_TITLE, strExtra);
//       bWriteVectorGISFile(VECTOR_PLOT_COAST, &VECTOR_PLOT_COAST_TITLE, strExtra);
//    }
// }
// #endif
