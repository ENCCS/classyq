#include "Sphere.hpp"

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace classyq {
Sphere::Sphere(double max_w, double r, Eigen::Vector3d c)
    : radius_(r), center_(c) {
  auto area = surface_of_sphere(radius_);
  // number of points in Leopardi partitioning of the unit sphere.
  N_ = static_cast<size_t>(std::ceil(area / max_w));

  // Leopardi partition of the unit sphere at the origin.
  auto [w_0, points_sph] = leopardi_partition(N_);

  // set w_0_ to the value computed
  w_0_ = w_0;

  // compute normal vectors at the EQ points by transforming EQ points from
  // spherical to Cartesian coordinates on the unit sphere
  normals_ = spherical_to_cartesian(points_sph);

  // from the normal vectors at EQ points (EQ points on the unit sphere)
  // we compute EQ points for this sphere as an affine transformation
  Eigen::Transform<double, 3, Eigen::Affine> T =
      Eigen::Translation3d(center_) * Eigen::Scaling(radius_);
  points_ = T * normals_;

  // compute the self-factors from the EQ points on the unit sphere
  self_potentials_ = Eigen::VectorXd(N_).setConstant(4.0 * M_PI);
  self_fields_ = Eigen::VectorXd(N_).setConstant(-2.0 * M_PI);

  // compute self-potential and self-field factors for the unit sphere
  for (auto i = 0; i < N_; ++i) {
    auto p_i = normals_.col(i);
    for (auto j = 0; j < N_; ++j) {
      auto p_j = normals_.col(j);
      if (i != j) {
        auto dist = (p_i - p_j).norm();
        self_potentials_(i) -= w_0_ / dist;
        self_fields_(i) -= (w_0_ * (p_i - p_j).dot(p_j)) / std::pow(dist, 3);
      }
    }
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
