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

  /** Return number of quadrature points. */
  auto size() const -> size_t { return N_; }

  /** Return number of spheres. */
  auto n_spheres() const -> size_t { return spheres_.size(); }
  /** Return list of spheres. */
  auto spheres() const -> std::vector<Sphere> { return spheres_; }

  /** Return a sphere in list.
   * @param[in] A index of sphere.
   */
  auto spheres(size_t A) const -> Sphere { return spheres_[A]; }

  /** Return threshold for weight of point. */
  auto threshold() const -> double { return threshold_; }

  /** Return vector of quadrature weights. */
  auto weights() const -> Eigen::VectorXd { return weights_; }

  /** Return list with number of points in quadrature belonging to each sphere. */
  auto points_per_sphere() const -> std::vector<size_t> { return points_per_sphere_; }

  /** Return number of points in quadrature belonging to given sphere.
   * @param[in] A index of sphere.
   */
  auto points_per_sphere(size_t A) const -> size_t { return points_per_sphere_[A]; }

  /** Return quadrature points. */
  auto points() const -> Eigen::Matrix3Xd { return points_; }

  /** Return given quadrature points.
   * @param[in] i index of point.
   */
  auto points(size_t i) const -> Eigen::Vector3d { return points_.col(i); }

  /** Return normal vectors at quadrature points. */
  auto normals() const -> Eigen::Matrix3Xd { return normals_; }

  /** Return given normal vector at quadrature points.
   * @param[in] i index of point.
   */
  auto normals(size_t i) const -> Eigen::Vector3d { return normals_.col(i); }

  /** Return surface area of cavity. */
  auto area() const -> double { return weights_.sum(); }

  /** Return surface area of cavity. */
  auto volume() const -> double;

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
