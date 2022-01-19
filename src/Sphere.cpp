#include "Sphere.hpp"

#include <cmath>

#include <Eigen/Core>

namespace classyq {
Sphere::Sphere(double max_w, double r, Eigen::Vector3d c)
    : radius_(r), center_(c) {
  auto area = surface_of_sphere(radius_);
  // number of points in Leopardi partitioning of the unit sphere.
  N_ = static_cast<size_t>(std::ceil(area / max_w));
  // adjust w_0_ so that N * w_0_ == area
  w_0_ = area / N_;

  // Leopardi partition of the unit sphere.
  auto [thetas, phis] = leopardi_partition(N_);
  // transfrom points from spherical to Cartesian coordinates
  points_ = Eigen::Matrix3Xd::Zero(3, N_);
  for (auto i = 0; i < N_; ++i) {
    points_.col(i) =
        spherical_to_cartesian(thetas(i), phis(i), radius_, center_);
  }
}

auto Sphere::switching(const Eigen::Vector3d &p, double rho) const -> double {
  // signed normal distance
  auto s = (p - center_).norm() - radius_;
  // penetration distance
  auto x = s / rho;

  return 0.5 * (1.0 + std::erf(x));
}
} // namespace classyq
