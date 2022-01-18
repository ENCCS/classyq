#pragma once

#include <tuple>

#include <Eigen/Core>

/** Compute EQ partitioning of the unit sphere.
 *
 * @param[in] N number of regions/points in the partition.
 * @return The EQ partition in spherical coordinates as tuple of vectors of
 * polar (\f$\theta\f$) and azimuthal (\f$\phi\f$) angles.
 *
 * The EQ partition subdivides the unit sphere into N regions of equal area.
 * This is a re-implementation (from Matlab) of the algorithm first described by
 * Leopardi in:
 *
 * Leopardi, P. A Partition of the Unit Sphere into Regions of Equal Area and
 * Small Diameter. In Electronic Transactions on Numerical Analysis.
 *
 * The original Matlab code is distributed under the MIT license.
 */
auto leopardi_partition(size_t N)
    -> std::tuple<Eigen::VectorXd, Eigen::VectorXd>;
