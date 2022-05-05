#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>
#include <fmt/format.h>

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
 */
class TsLess final
{
private:
  /** Total number of quadrature points. */
  size_t N_{ 0 };

  /** Quadrature weight threshold. */
  double threshold_{ std::numeric_limits<double>::epsilon() };

  /** List of spheres */
  std::vector<Sphere> spheres_;

  /** Weight of quadrature points. */
  Eigen::VectorXd weights_;

  /** Quadrature points. */
  Eigen::Matrix3Xd points_;

  /** Normal vectors at quadrature points. */
  Eigen::Matrix3Xd normals_;

  /** List of points per sphere in the final quadrature. */
  std::vector<size_t> points_per_sphere_;

  /** List of exposed areas per sphere. */
  std::vector<double> exposed_area_;

public:
  TsLess(const std::vector<Sphere>& ss, double threshold);

  auto size() const -> size_t { return N_; }

  auto threshold() const -> double { return threshold_; }

  auto weights() const -> Eigen::VectorXd { return weights_; }

  auto points() const -> Eigen::Matrix3Xd { return points_; }

  auto points(size_t i) const -> Eigen::Vector3d { return points_.col(i); }

  auto normals() const -> Eigen::Matrix3Xd { return normals_; }

  auto normals(size_t i) const -> Eigen::Vector3d { return normals_.col(i); }

  auto area() const -> double { return weights_.sum(); }

  /** String representation of object. */
  auto str() const -> std::string;
};
} // namespace classyq

template<>
struct fmt::formatter<classyq::TsLess> : fmt::formatter<std::string_view>
{
  template<typename FormatContext>
  auto format(const classyq::TsLess& obj, FormatContext& ctx) -> decltype(ctx.out())
  {
    return fmt::formatter<std::string_view>::format(obj.str(), ctx);
  }
};
