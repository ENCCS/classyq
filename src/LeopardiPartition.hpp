#pragma once

#include <tuple>

#include <Eigen/Core>

namespace classyq {
/** Compute EQ partitioning of the unit sphere.
 *
 * @param[in] N number of regions/points in the partition.
 * @return The weight and points of the EQ partition for the unit sphere at the origin.
 * The points are given in spherical coordinates as (2 x N) matrix of
 * polar (\f$\theta\f$) and azimuthal (\f$\phi\f$) angles.
 *
 * The EQ partition subdivides the unit sphere into N regions of equal area.
 * This is a re-implementation (from Matlab) of the algorithm first described by
 * Leopardi in \cite Leopardi2006-lh
 *
 * The original Matlab code is distributed under the MIT license.
 */
auto leopardi_partition(size_t N) -> std::tuple<double, Eigen::Matrix2Xd>;
} // namespace classyq
