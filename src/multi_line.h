/*!
   \class CGeomMultiLine
   \brief Geometry class used to represent co-incident lines (for profiles/polygon-to-polygon boundaries)
   \details TODO 001 This is a more detailed description of the CGeomMultiLine class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file multi_line.h
   \brief Contains CGeomMultiLine definitions
*/

#ifndef MULTILINE_H
#define MULTILINE_H
/* ===============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
#include <vector>
using std::vector;

#include <utility>
using std::pair;
using std::make_pair;

#include "2d_point.h"
#include "line.h"

class CGeomMultiLine : public CGeomLine
{
 private:
   //! A vector of line segments, each element is a vector of pairs. The first of the pair is a co-incident profile number, the second is that profile's 'own' line segment number
   vector<vector<pair<int, int>>> m_prVVLineSegment;

 protected:
 public:
   CGeomMultiLine(void);
   ~CGeomMultiLine(void) override;

   vector<CGeom2DPoint>& pGetPoints(void);
   // void SetPoints(vector<CGeom2DPoint> const&);

   void AppendLineSegment(void);
   void AppendLineSegment(vector<pair<int, int>>*);
   // void AppendLineSegmentAndInherit(void);
   int nGetNumLineSegments(void) const;
   void TruncateLineSegments(int const);
   void InsertLineSegment(int const);
   vector<vector<pair<int, int>>> prVVGetAllLineSegAfter(int const);
   // void RemoveLineSegment(int const);

   void AppendCoincidentProfileToLineSegments(pair<int, int> const);
   void AddCoincidentProfileToExistingLineSegment(int const, int const, int const);
   vector<pair<int, int>>* pprVGetPairedCoincidentProfilesForLineSegment(int const);
   int nGetCoincidentProfileForLineSegment(int const, int const) const;
   int nGetNumCoincidentProfilesInLineSegment(int const);
   bool bFindProfileInCoincidentProfilesOfLastLineSegment(int const);
   // bool bFindProfileInCoincidentProfilesOfLineSegment(int const, int const);
   bool bFindProfileInCoincidentProfiles(int const);
   void GetMostCoastwardSharedLineSegment(int const, int&, int &);

   int nGetProf(int const, int const) const;
   int nGetProfsLineSeg(int const, int const) const;
   void SetProfsLineSeg(int const, int const, int const);

   // int nFindProfilesLastSeg(int const) const;
};
#endif // MULTILINE_H
