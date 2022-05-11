/*
 * ClassyQ, classical polarizable solvent models.
 * Copyright (C) 2022 Roberto Di Remigio Eikås and contributors.
 *
 * This file is part of ClassyQ.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * ClassyQ library, see: <https://classyq.readthedocs.io/>
 */

#pragma once

#include <string>
#include <vector>

#include <fmt/format.h>

#include <Eigen/Core>

namespace classyq {
/**
 */
class Sphere final
{
private:
  /** Whether the sphere is centered on an atom.
   *
   * Spheres not attached to an atom do not contribute to the molecular gradient.
   */
  bool atom_centered_{ true };

  /** Number of EQ points. */
  size_t N_{ 0 };

  /** Weight of EQ points **for the unit sphere**.
   *
   * This is set to \f$ \omega = \frac{4\pi}{N} \f$, while the weight is \f$w =
   * \omega R^{2}\f$. The former is used in the definition of the self-potential
   * and self-field factors, while the latter is the weight used in the
   * quadrature of the interlocking-sphere cavity.
   */
  double omega_{ 0.0 };

  /** Radius. */
  double radius_{ 0.0 };

  /** Center. */
  Eigen::Vector3d center_{ 0.0, 0.0, 0.0 };

  /** EQ points, in Cartesian coordinates, **for the unit sphere**.
   *
   * @note these are equivalent to the normal vectors at the points on the
   * non-unit sphere.
   */
  Eigen::Matrix3Xd points_;

public:
  /** Constructor from maximum quadrature weight.
   *
   * @param[in] max_w maximum quadrature weight (inverse of the density of points)
   * @param[in] r radius of the sphere.
   * @param[in] c center of the sphere.
   * @param[in] on_atom whether the sphere is centered on an atom.
   * Spheres not centered on atoms do not contribute to the molecular gradient.
   *
   * TODO convert it to density of points (G16 uses 5 points per Ang^2)
   */
  Sphere(double max_w, double r, const Eigen::Vector3d& c, bool on_atom = true);

  /** Switching function at given point.
   *
   * @param[in] s evaluation point.
   * @param[in] rho interaction radius of point.
   *
   * The switching function is a smoothed Heaviside step function:
   *
   * - 0 if point \f$j\f$ on sphere \f$\alpha\f$ is inside sphere \f$\beta\f$ (this sphere),
   * - 1 if it's outside.
   *
   * Formally:
   *
   * \f[
   *   \vartheta(d_{\beta}^{\alpha j})
   *   =
   *   \frac{1}{2}\left[ 1 + \erf\left(d_{\beta}^{\alpha j}\right) \right],
   * \f]
   *
   * where \f$d_{\beta}^{\alpha j}\f$ is the *penetration distance*:
   *
   * \f[
   *   d_{\beta}^{\alpha j}
   *   =
   *   \frac{|\mathbf{s}_{\alpha j} - \mathbf{c}_{\beta}| - R_{\beta}}{\rho},
   * \f]
   *
   * defined in terms of the *signed normal distance* and the *interaction
   * radius* \f$\rho\f$. The latter is a function of the given evaluation point.
   */
  auto switching(const Eigen::Vector3d& s, double rho) const -> double;

  /** Center of sphere. */
  auto center() const -> Eigen::Vector3d { return center_; }

  /** Radius of sphere. */
  auto radius() const -> double { return radius_; }

  /** Number of quadrature points. */
  auto size() const -> size_t { return N_; }

  /** Weight of EQ points. */
  auto weight() const -> double { return omega_ * std::pow(radius_, 2); }

  /** Weight of EQ points **on the unit sphere** */
  auto omega() const -> double { return omega_; }

  /** EQ points. */
  auto points() const -> Eigen::Matrix3Xd { return points_; }

  auto is_atom_centered() const -> bool { return atom_centered_; }

  /** String representation of object. */
  auto str() const -> std::string;
};

using NeighborsList = std::vector<std::vector<size_t>>;

/** Computes neighbor list of a set of spheres.
 *
 * \param[in] ss the list of spheres.
 * \return the neighbor list as a list of lists.
 */
auto
neighbors_list(const std::vector<Sphere>& ss) -> NeighborsList;
} // namespace classyq

template<>
struct fmt::formatter<classyq::Sphere> : fmt::formatter<std::string_view>
{
  template<typename FormatContext>
  auto format(const classyq::Sphere& obj, FormatContext& ctx) -> decltype(ctx.out())
  {
    return fmt::formatter<std::string_view>::format(obj.str(), ctx);
  }
};
