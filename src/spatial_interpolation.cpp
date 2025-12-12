/*!
   \file spatial_interpolation.cpp
   \brief Spatial Interpolation Using k-Nearest Neighbors and Inverse Distance Weighting
   \details TODO
   \author Wilf Chun
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/

//===============================================================================================================================
//! Spatial Interpolation Using k-Nearest Neighbors and Inverse Distance Weighting
//!
//! This file implements fast spatial interpolation using:
//! 1. k-d tree for efficient nearest neighbor search (O(log n) per query)
//! 2. Inverse Distance Weighting (IDW) for smooth interpolation
//! 3. OpenMP parallelization for batch operations
//!
//! ALGORITHM OVERVIEW:
//! ------------------
//! For each query point (grid cell):
//!   1. Find k nearest input points using k-d tree
//!   2. Calculate weights based on inverse distance: w_i = 1 / dist_i^power
//!   3. Interpolated value = Σ(w_i * value_i) / Σ(w_i)
//!
//! KEY PARAMETERS TO TUNE:
//! ----------------------
//! - k_neighbors: Number of nearest points to use (typically 8-20)
//!   * Higher = smoother, more averaging
//!   * Lower = follows local variations more closely
//!
//! - power: Exponent for inverse distance weighting (typically 1.0-4.0)
//!   * Higher = nearby points have more influence (sharper transitions)
//!   * Lower = distant points have more influence (smoother transitions)
//!   * power=2.0 is optimized for performance (avoids pow() function)
//!
//! PERFORMANCE OPTIMIZATIONS:
//! -------------------------
//! - Special fast path for power=2.0 (uses squared distance directly, avoids sqrt/pow)
//! - OpenMP parallelization with guided scheduling
//! - Thread-local buffers to avoid allocation overhead
//! - Shared k-d tree between X and Y interpolation (DualSpatialInterpolator)
//!
//! USAGE EXAMPLE:
//! --------------
//! std::vector<Point2D> input_points = {{0,0}, {10,0}, {5,10}};
//! std::vector<double> input_values = {1.0, 2.0, 1.5};
//! SpatialInterpolator interp(input_points, input_values, 12, 2.0);
//! double result = interp.Interpolate(5.0, 5.0);  // Interpolate at (5,5)
//!
//===============================================================================================================================
#include <stdexcept>

#include "spatial_interpolation.h"
#include "cme.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//===============================================================================================================================
//! Constructor: Build interpolator from points and values
//!
//! @param points      Input point coordinates (x, y)
//! @param values      Values at those points
//! @param k_neighbors Number of nearest neighbors to use (default: 12)
//! @param power       IDW power parameter (default: 2.0)
//===============================================================================================================================
SpatialInterpolator::SpatialInterpolator(std::vector<Point2D> const& points,
                                         std::vector<double> const& values,
                                         int k_neighbors,
                                         double power)
   : m_values(values), m_kdtree(nullptr), m_k_neighbors(k_neighbors), m_power(power), m_owns_kdtree(true)
{
   if (points.size() != values.size())
      throw std::invalid_argument("Points and values must have same size");

   if (points.empty())
      throw std::invalid_argument("Cannot create interpolator with empty data");

   // Copy input points into point cloud structure
   m_cloud.pts = points;

   // Build k-d tree for fast spatial queries
   // max_leaf_size=10 is a good balance between build time and query speed
   m_kdtree = new KDTree(2, m_cloud, {10 /* max leaf size */});
   m_kdtree->buildIndex();
}

SpatialInterpolator::SpatialInterpolator(PointCloud const& cloud,
                                         KDTree* kdtree,
                                         std::vector<double> const& values,
                                         int k_neighbors,
                                         double power)
   : m_cloud(cloud), m_values(values), m_kdtree(kdtree), m_k_neighbors(k_neighbors), m_power(power), m_owns_kdtree(false)
{
}

SpatialInterpolator::~SpatialInterpolator()
{
   if (m_owns_kdtree)
      delete m_kdtree;
}

//===============================================================================================================================
//! Interpolate at a single query point
//!
//! ALGORITHM:
//! 1. Find k nearest neighbors using k-d tree
//! 2. Calculate weight for each: w_i = 1 / distance_i^power
//! 3. Return weighted average: Σ(w_i * value_i) / Σ(w_i)
//!
//! SPECIAL CASES:
//! - If query point coincides with input point (dist < EPSILON): return exact value
//! - If power=2.0: optimized path using squared distances (faster)
//!
//! @param x  X coordinate of query point
//! @param y  Y coordinate of query point
//! @return Interpolated value
//===============================================================================================================================
double SpatialInterpolator::Interpolate(double x, double y) const
{
   // Prepare query point
   double const query_pt[2] = {x, y};

   // Find k nearest neighbors (or all points if fewer than k exist)
   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());
   std::vector<unsigned int> indices(k);
   std::vector<double> sq_dists(k);  // Squared distances (faster than actual distances)

   long unsigned int num_found = m_kdtree->knnSearch(query_pt, k,
                                                 indices.data(),
                                                 sq_dists.data());

   if (num_found == 0)
      throw std::runtime_error("knnSearch found no neighbors");

   // SPECIAL CASE: Query point coincides with an input point
   // Return exact value to avoid division by zero
   if (sq_dists[0] < EPSILON)
      return m_values[indices[0]];

   // Inverse Distance Weighting (IDW)
   // Formula: result = Σ(w_i * v_i) / Σ(w_i) where w_i = 1/dist_i^power
   double sum_weights = 0.0;
   double sum_weighted_values = 0.0;

   if (bFPIsEqual(m_power, 2.0, TOLERANCE))
   {
      // *** OPTIMIZED PATH for power=2.0 ***
      // Since weight = 1/dist^2 and we have sq_dist = dist^2,
      // we can use: weight = 1/sq_dist
      // This avoids both sqrt() and pow() calls
      for (size_t i = 0; i < num_found; i++)
      {
         double weight = 1.0 / sq_dists[i];  // 1/dist^2 = 1/sq_dist
         sum_weights += weight;
         sum_weighted_values += weight * m_values[indices[i]];
      }
   }
   else
   {
      // *** GENERAL CASE for arbitrary power ***
      // Need to calculate actual distance and apply pow()
      for (size_t i = 0; i < num_found; i++)
      {
         double dist = std::sqrt(sq_dists[i]);
         double weight = 1.0 / std::pow(dist, m_power);
         sum_weights += weight;
         sum_weighted_values += weight * m_values[indices[i]];
      }
   }

   return sum_weighted_values / sum_weights;
}

void SpatialInterpolator::Interpolate(std::vector<Point2D> const& query_points,
                                      std::vector<double>& results) const
{
   results.resize(query_points.size());
   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());

#ifdef _OPENMP
   #pragma omp parallel
   {
      // Thread-local storage to avoid allocation overhead
      std::vector<unsigned int> indices(k);
      std::vector<double> sq_dists(k);

      #pragma omp for schedule(guided, 128) nowait
      for (size_t i = 0; i < query_points.size(); i++)
      {
         double const query_pt[2] = {query_points[i].x, query_points[i].y};

         long unsigned int num_found = m_kdtree->knnSearch(query_pt, k,
                                                       indices.data(),
                                                       sq_dists.data());

         if (num_found == 0)
         {
            results[i] = 0.0;  // or m_dMissingValue
            continue;
         }

         // Check if we're exactly on a data point
         if (sq_dists[0] < EPSILON)
         {
            results[i] = m_values[indices[0]];
            continue;
         }

         // IDW calculation - optimized for power=2.0
         double sum_weights = 0.0;
         double sum_weighted_values = 0.0;

         if (bFPIsEqual(m_power, 2.0, TOLERANCE))
         {
            for (size_t j = 0; j < num_found; j++)
            {
               double weight = 1.0 / sq_dists[j];
               sum_weights += weight;
               sum_weighted_values += weight * m_values[indices[j]];
            }
         }
         else
         {
            for (size_t j = 0; j < num_found; j++)
            {
               double dist = std::sqrt(sq_dists[j]);
               double weight = 1.0 / std::pow(dist, m_power);
               sum_weights += weight;
               sum_weighted_values += weight * m_values[indices[j]];
            }
         }

         results[i] = sum_weighted_values / sum_weights;
      }
   }
#else
   // Serial fallback - still optimized with reused buffers
   std::vector<unsigned int> indices(k);
   std::vector<double> sq_dists(k);

   for (size_t i = 0; i < query_points.size(); i++)
   {
      results[i] = Interpolate(query_points[i].x, query_points[i].y);
   }
#endif
}

//===============================================================================================================================
//! DualSpatialInterpolator: Optimized interpolation for paired X/Y values
//!
//! This class is optimized for interpolating wave properties that have both X and Y components
//! (e.g., wave height X and wave height Y). It's more efficient than using two separate
//! interpolators because:
//! 1. Builds only ONE k-d tree (shared for both X and Y)
//! 2. Does ONE nearest neighbor search per query point (not two)
//! 3. Interpolates both components simultaneously
//!
//! USAGE:
//! ------
//! DualSpatialInterpolator interp(points, values_x, values_y, 12, 2.0);
//! interp.Interpolate(query_points, results_x, results_y);
//!
//===============================================================================================================================

//===============================================================================================================================
//! Constructor: Build dual interpolator for X and Y values
//!
//! @param points      Input point coordinates (x, y)
//! @param values_x    X-component values at those points
//! @param values_y    Y-component values at those points
//! @param k_neighbors Number of nearest neighbors to use (default: 12)
//! @param power       IDW power parameter (default: 2.0)
//===============================================================================================================================
DualSpatialInterpolator::DualSpatialInterpolator(std::vector<Point2D> const& points,
                                                 std::vector<double> const& values_x,
                                                 std::vector<double> const& values_y,
                                                 int k_neighbors,
                                                 double power)
   : m_values_x(values_x), m_values_y(values_y), m_kdtree(nullptr), m_k_neighbors(k_neighbors), m_power(power)
{
   if (points.size() != values_x.size() || points.size() != values_y.size())
      throw std::invalid_argument("Points and values must have same size");

   if (points.empty())
      throw std::invalid_argument("Cannot create interpolator with empty data");

   // Copy input points
   m_cloud.pts = points;

   // Build k-d tree (shared for both X and Y interpolation)
   m_kdtree = new KDTree(2, m_cloud, {10 /* max leaf size */});
   m_kdtree->buildIndex();
}

DualSpatialInterpolator::~DualSpatialInterpolator()
{
   delete m_kdtree;
}

void DualSpatialInterpolator::InterpolatePoint(double x, double y,
                                               double& result_x, double& result_y,
                                               std::vector<unsigned int>& indices,
                                               std::vector<double>& sq_dists) const
{
   double const query_pt[2] = {x, y};
   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());

   long unsigned int num_found = m_kdtree->knnSearch(query_pt, k,
                                                 indices.data(),
                                                 sq_dists.data());

   if (num_found == 0)
   {
      result_x = result_y = 0.0;
      return;
   }

   // Check if we're exactly on a data point
   if (sq_dists[0] < EPSILON)
   {
      result_x = m_values_x[indices[0]];
      result_y = m_values_y[indices[0]];
      return;
   }

   // IDW for both X and Y simultaneously - optimized for power=2.0
   double sum_weights = 0.0;
   double sum_weighted_x = 0.0;
   double sum_weighted_y = 0.0;

   if (bFPIsEqual(m_power, 2.0, TOLERANCE))
   {
      for (size_t i = 0; i < num_found; i++)
      {
         double weight = 1.0 / sq_dists[i];
         sum_weights += weight;
         sum_weighted_x += weight * m_values_x[indices[i]];
         sum_weighted_y += weight * m_values_y[indices[i]];
      }
   }
   else
   {
      for (size_t i = 0; i < num_found; i++)
      {
         double dist = std::sqrt(sq_dists[i]);
         double weight = 1.0 / std::pow(dist, m_power);
         sum_weights += weight;
         sum_weighted_x += weight * m_values_x[indices[i]];
         sum_weighted_y += weight * m_values_y[indices[i]];
      }
   }

   result_x = sum_weighted_x / sum_weights;
   result_y = sum_weighted_y / sum_weights;
}

void DualSpatialInterpolator::Interpolate(std::vector<Point2D> const& query_points,
                                          std::vector<double>& results_x,
                                          std::vector<double>& results_y) const
{
   size_t const n = query_points.size();
   results_x.resize(n);
   results_y.resize(n);

   // Early exit for empty query
   if (n == 0) return;

   size_t const k = std::min((size_t) m_k_neighbors, m_cloud.pts.size());

#ifdef _OPENMP
   // Use guided scheduling to reduce overhead
   #pragma omp parallel
   {
      // Thread-local storage
      std::vector<unsigned int> indices(k);
      std::vector<double> sq_dists(k);

      #pragma omp for schedule(guided, 128) nowait
      for (size_t i = 0; i < query_points.size(); i++)
      {
         InterpolatePoint(query_points[i].x, query_points[i].y,
                         results_x[i], results_y[i],
                         indices, sq_dists);
      }
   }
#else
   // Serial fallback
   std::vector<unsigned int> indices(k);
   std::vector<double> sq_dists(k);

   for (size_t i = 0; i < query_points.size(); i++)
   {
      InterpolatePoint(query_points[i].x, query_points[i].y,
                      results_x[i], results_y[i],
                      indices, sq_dists);
   }
#endif
}
