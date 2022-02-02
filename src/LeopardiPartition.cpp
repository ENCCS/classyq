#include "LeopardiPartition.hpp"

#include <cmath>
#include <numeric>
#include <tuple>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <Eigen/Core>

namespace classyq {
auto leopardi_partition(size_t N) -> std::tuple<double, Eigen::Matrix2Xd> {
  // IDEAL area of each region
  auto V_R = (4.0 * M_PI) / N;
  SPDLOG_TRACE("ideal area of regions: V_R = {}", V_R);

  // colatitude of the north pole spherical cap (spherical radius of spherical
  // cap of given area)
  auto theta_c = 2.0 * std::asin(std::sqrt(1.0 / N));
  // the colatitude of the North pole is theta_c, the colatitude of the South
  // pole is M_PI - theta_c

  // IDEAL collar angle
  auto delta_I = std::sqrt(V_R);

  // IDEAL number of collars between polar caps
  // this is n_I in the paper
  auto n_I = static_cast<size_t>(
      std::max(1.0, std::round((M_PI - 2 * theta_c) / delta_I)));
  SPDLOG_TRACE("ideal number of collars between polar caps: n_I = {}", n_I);

  // FITTING collar angle
  auto delta_F = (M_PI - 2.0 * theta_c) / n_I;

  // total number of regions. We always have at least 2 (the poles)
  auto n_regions = 2;

  // array holding number of regions in each collar
  // y[0] = 1 (only one region in the North pole)
  // y[n_I+1] = 1 (only one region in the South pole)
  auto y = std::vector<size_t>(n_I + 2, 1);
  // accumulated difference between ideal and actual number of regions.
  auto diff = 0.0;
  // fill the y array: for each of the n_I collars, compute the ideal number of
  // regions.  This is not integer, so we add the difference accumulated so far
  // to it and round to integer.  Accumulating the difference between actual and
  // ideal guarantees that we get the requested number of regions!
  for (auto i = 1; i <= n_I; ++i) {
    // the paper is generic in the number of dimensions of the sphere,
    // we specialize for the 2-sphere and use a different formula for the ideal
    // number of regions (saves some FLOPs).
    // TODO write formulas here
    auto ideal =
        N * std::sin(theta_c + (i - 0.5) * delta_F) * std::sin(delta_F / 2.0);
    y[i] = std::round(ideal + diff);
    diff += ideal - y[i];
    n_regions += y[i];
  }
  SPDLOG_TRACE("number of regions per collar: y = {}", y);
  if (n_regions != N) {
    SPDLOG_WARN("computed number of regions ({}) differs from requested ({})",
                n_regions, N);
  }

  // colatitudes of each zone
  auto theta = std::vector<double>(n_I + 2, 0.0);
  // first element is the colatitude of the north pole spherical cap
  theta.front() = theta_c;
  // one-before-last element is the colatitudes of the south pole spherical cap:
  // theta[n_I] = M_PI - theta_c
  // last element is the south pole itself
  theta.back() = M_PI;

  auto acc = static_cast<double>(y[0]);
  for (auto i = 1; i <= n_I; ++i) {
    acc += y[i];
    theta[i] = 2.0 * std::asin(std::sqrt(acc / N));
  }
  SPDLOG_TRACE("colatitudes of each zone: theta = {}", theta);

  // result array with polar (row 0) and azimuthal (row 1) angles
  Eigen::Matrix2Xd points(2, n_regions);
  // north pole
  points.col(0) = Eigen::Vector2d::Zero();
  // south pole
  points.col(n_regions - 1) << M_PI, 0.0;

  // loop over collars, excluding north and south poles
  auto start = y[0];
  auto offset = 0.0;
  for (auto i = 1; i <= n_I; ++i) {
    // a_top is the colatitude of the top of the current collar.
    auto a_top = theta[i - 1];
    // a_bot is the colatitude of the bottom of the current collar.
    auto a_bot = theta[i];
    auto n_regions = y[i];

    auto a_point = (a_top + a_bot) / 2.0;

    // all regions in this collar have the same polar angle
    points(0, Eigen::seqN(start, n_regions)).setConstant(a_point);

    // each region is a circle, which we can partition trivially
    auto inc = 2.0 * M_PI / n_regions;
    for (auto j = 1; j <= n_regions; ++j) {
      auto x = j * inc - M_PI / n_regions;
      // azimuthal angles
      // we "twist" the azimuthal partition such that points in different
      // regions do not align.
      points(1, j + (start - 1)) = std::fmod(x + 2 * M_PI * offset, 2 * M_PI);
    }
    // update the azimuthal "twist", as follows:
    //  1. half the difference between a "twist" of one sector on each of bottom
    //  and top. This brings the centre points into alignment.
    //  2. A rotation which will maximize the minimum angle between points on
    //  the two circles.
    // We accumulate this offset, and force it to be a number between 0 and 1.
    auto n_regions_next = y[i + 1];
    offset += 0.5 * (1.0 / n_regions_next - 1.0 / n_regions) +
              std::gcd(n_regions, n_regions_next) /
                  (2.0 * n_regions * n_regions_next);
    offset -= std::floor(offset);

    start += n_regions;
  }

  return {V_R, points};
}
} // namespace classyq
