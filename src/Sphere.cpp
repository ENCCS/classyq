#include "Sphere.hpp"

#include <cmath>
#include <tuple>

#include <Eigen/Core>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include "LeopardiPartition.hpp"

namespace classyq {
Sphere::Sphere(double max_w, double r, Eigen::Vector3d c)
    : radius_(r), center_(c) {
  auto area = surface_of_sphere(radius_);
  // number of points in Leopardi partitioning of the unit sphere.
  N_ = static_cast<size_t>(std::ceil(area / max_w));

  // Leopardi partition of the unit sphere at the origin.
  Eigen::Matrix2Xd points_sph;
  std::tie(w_0_, points_sph) = leopardi_partition(N_);

  // compute normal vectors at the EQ points by transforming EQ points from
  // spherical to Cartesian coordinates on the unit sphere
  points_ = spherical_to_cartesian(points_sph);

  // compute the self-factors from the EQ points on the unit sphere
  self_potentials_ = Eigen::VectorXd(N_).setConstant(4.0 * M_PI);
  self_fields_ = Eigen::VectorXd(N_).setConstant(-2.0 * M_PI);

  // compute self-potential and self-field factors for the unit sphere
  for (auto i = 0; i < N_; ++i) {
    auto p_i = points_.col(i);
    for (auto j = 0; j < N_; ++j) {
      auto p_j = points_.col(j);
      if (i != j) {
        auto dist = (p_i - p_j).norm();
        self_potentials_(i) -= w_0_ / dist;
        self_fields_(i) -= (w_0_ * (p_i - p_j).dot(p_j)) / std::pow(dist, 3);
      }
    }
  }
}
} // namespace classyq
