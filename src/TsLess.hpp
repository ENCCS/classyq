#pragma once

#include <numeric>
#include <tuple>
#include <vector>

#include <Eigen/Core>

#include "Sphere.hpp"

namespace classyq {
/**
 *
 * This is a modified implementation of the cavity partition first described in
 * \cite Pomelli2004-lb
 * At variance with the work of Pomelli:
 *
 * 1. We define the sphere-sphere switching using the error function, instead of
 * a polynomial.
 * 2. The diagonal factors for the collocation of the boundary integral
 * operators are defined as in \cite Scalmani2010-tw.
 * 3. Finally, the geometric derivatives of the cavity are computed with
 * automatic differentiation, instead of being coded explicitly.
 */
class TsLess final {
private:
  size_t N_{0};
  Eigen::VectorXd weights_;
  Eigen::Matrix3Xd normals_;
  Eigen::Matrix3Xd points_;
  Eigen::VectorXd self_potentials_;
  Eigen::VectorXd self_fields_;

  std::vector<size_t> points_per_sphere_;
  std::vector<double> weights_0_;
  std::vector<double> radii_;

public:
  TsLess(const std::vector<Sphere> &ss,
         double threshold = std::numeric_limits<double>::epsilon());

  auto size() const -> size_t { return N_; }

  auto weights() const -> Eigen::VectorXd { return weights_; }

  auto points() const -> Eigen::Matrix3Xd { return points_; }
  auto points(size_t i) const -> Eigen::Vector3d { return points_.col(i); }

  auto normals() const -> Eigen::Matrix3Xd { return normals_; }
  auto normals(size_t i) const -> Eigen::Vector3d { return normals_.col(i); }

  auto area() const -> double { return weights_.sum(); }

  auto fs() const -> Eigen::VectorXd { return self_potentials_; }
  auto gs() const -> Eigen::VectorXd { return self_fields_; }
  auto self_factors() const -> std::tuple<Eigen::VectorXd, Eigen::VectorXd>;

  /**
   * @{ Iterators
   */
  /**@}*/
};

using NeighborsList = std::vector<std::vector<size_t>>;

auto neighbors_list(const std::vector<Sphere> &spheres) -> NeighborsList;
} // namespace classyq
