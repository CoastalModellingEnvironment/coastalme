/*!
   \file spatial_interpolation.ch
   \brief Spatial Interpolation Using k-Nearest Neighbors and Inverse Distance Weighting
   \details TODO
   \autho Wilf Chun
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
*/


#ifndef SPATIAL_INTERPOLATION_H
#define SPATIAL_INTERPOLATION_H

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/

//===============================================================================================================================
//! Spatial Interpolation Using k-Nearest Neighbors and Inverse Distance Weighting
//!
//! This header defines classes for fast spatial interpolation using:
//! - k-d tree spatial indexing (nanoflann library)
//! - Inverse Distance Weighting (IDW)
//! - OpenMP parallelization
//!
//! CLASSES:
//! --------
//! - Point2D: Simple 2D point structure
//! - PointCloud: k-d tree adaptor for nanoflann
//! - SpatialInterpolator: Single-value interpolation
//! - DualSpatialInterpolator: Optimized paired X/Y interpolation (recommended for wave data)
//!
//! QUICK START:
//! ------------
//! vector<Point2D> points = {{0,0}, {10,0}, {5,10}};
//! vector<double> values_x = {1.0, 2.0, 1.5};
//! vector<double> values_y = {0.5, 1.0, 0.8};
//! DualSpatialInterpolator interp(points, values_x, values_y, 12, 2.0);
//! vector<Point2D> query = {{5,5}};
//! vector<double> result_x, result_y;
//! interp.Interpolate(query, result_x, result_y);
//!
//! TUNING PARAMETERS:
//! -----------------
//! - k_neighbors (default 12): Number of nearest points to use
//!   * Increase for smoother results
//!   * Decrease for more local detail
//! - power (default 2.0): IDW power exponent
//!   * Increase for sharper transitions (nearby points dominate)
//!   * Decrease for smoother transitions (distant points have more influence)
//!
//===============================================================================================================================
#include <vector>
using std::vector;

#include <cmath>
#include <algorithm>

#include "nanoflann.hpp"

//===============================================================================================================================
//! Simple 2D point structure
//===============================================================================================================================
struct Point2D
{
   double x, y;
   explicit Point2D(double x_ = 0, double y_ = 0) : x(x_), y(y_)
   {
   }
};

// Point cloud adaptor for nanoflann
struct PointCloud
{
   std::vector<Point2D> pts;

   inline size_t kdtree_get_point_count() const
   {
      return pts.size();
   }

   inline double kdtree_get_pt(size_t const idx, size_t const dim) const
   {
      return (dim == 0) ? pts[idx].x : pts[idx].y;
   }

   template<class BBOX>
   bool kdtree_get_bbox(BBOX&) const
   {
      return false;
   }
};

// Type definitions for nanoflann
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>;

// Main interpolator class
class SpatialInterpolator
{
   public:
   SpatialInterpolator(std::vector<Point2D> const& points, std::vector<double> const& values, int k_neighbors = 12, double power = 2.0);

   ~SpatialInterpolator();

   // Interpolate at a single point
   double Interpolate(double x, double y) const;

   // Interpolate at multiple points (more efficient)
   void Interpolate(vector<Point2D> const& query_points, vector<double>& results) const;

   // Get the k-d tree (for sharing between interpolators)
   KDTree const* GetKDTree() const { return m_kdtree; }

   // Get the point cloud (for sharing between interpolators)
   PointCloud const& GetPointCloud() const { return m_cloud; }

   private:
   PointCloud m_cloud;
   vector<double> m_values;
   KDTree* m_kdtree;
   int m_k_neighbors;
   double m_power;
   bool m_owns_kdtree;

   static constexpr double EPSILON = 1e-10;

   // Private constructor for sharing k-d tree
   friend class DualSpatialInterpolator;
   SpatialInterpolator(PointCloud const& cloud, KDTree* kdtree, vector<double> const& values, int k_neighbors, double power);
};

// Optimized dual interpolator for X and Y values sharing same spatial points
class DualSpatialInterpolator
{
public:
   DualSpatialInterpolator(vector<Point2D> const& points, vector<double> const& values_x, vector<double> const& values_y, int k_neighbors = 12, double power = 2.0);

   ~DualSpatialInterpolator();

   // Interpolate both X and Y at multiple points (parallel optimized)
   void Interpolate(vector<Point2D> const& query_points, vector<double>& results_x, vector<double>& results_y) const;

private:
   PointCloud m_cloud;
   vector<double> m_values_x;
   vector<double> m_values_y;
   KDTree* m_kdtree;
   int m_k_neighbors;
   double m_power;

   static constexpr double EPSILON = 1e-10;

   // Helper for single point interpolation
   void InterpolatePoint(double x, double y, double& result_x, double& result_y, vector<unsigned int>& indices, vector<double>& sq_dists) const;
};

#endif      // SPATIAL_INTERPOLATION_H
