/*!
   \class CGeomMultiLine
   \brief Geometry class used to represent co-incident lines (for profiles/polygon-to-polygon boundaries)
   \details TODO 001 This is a more detailed description of the CGeomMultiLine class.
   \author David Favis-Mortlock
   \author Andres Payo
   \author Wilf Chun
   \date 2026
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
   //! Each CGeomMultiLine is a vector of vectors: a vector of line segments, with each line segment being a vector of pairs. The first element of each pair is a co-incident profile number, the second element is that profile's 'own' line segment number. Since the CGeomMultiLine is descended from a CGeomLine, the CGeomMultiLine also has a vector of points (m_VPoints). Each point is the start or finish of a line segment. Thus each CGeomMultiLine should have N line segments and N+1 points
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
   void InsertLineSegmentWithInheritance(int const);
   vector<vector<pair<int, int>>> prVVGetAllLineSegmentsAfter(int const);
   // void RemoveLineSegment(int const);

   void AppendPairToFinalLineSegment(pair<int, int> const);
   void AddCoincidentPairToExistingLineSegmentIfNotAlready(int const, int const, int const);

   vector<pair<int, int>>* pprVGetCoincidentPairsForLineSegment(int const);
   int pprVGetFirstFromCoincidentPairForLineSegment(int const, int const) const;

   int nGetNumCoincidentPairsForLineSegment(int const);

   bool bFindFirstFromCoincidentPairsInLastLineSegment(int const);
   // bool bFindFirstInCoincidentPairsOfLineSegment(int const, int const);
   bool bFindFirstInCoincidentPairs(int const);
   void SearchForLowestNumberedCoincidentLineSegments(int const, int&, int&);

   int nGetPairFirst(int const, int const) const;
   int nGetPairSecond(int const, int const) const;
   void SetPairSecond(int const, int const, int const);

   // int nFindLastSegForCoincPairFirst(int const) const;
};
#endif // MULTILINE_H
