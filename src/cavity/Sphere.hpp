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
   * This is set to \f$ w_{0} = \frac{4\pi}{N} \f$, while the weight is \f$w =
   * w_{0}R^{2}\f$. The former is used in the definition of the self-potential
   * and self-field factors, while the latter is the weight used in the
   * quadrature of the interlocking-sphere cavity.
   */
  double w_0_{ 0.0 };

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

  /** Self-potential factors **for the unit sphere**.
   *
   * The formula is given by Scalmani and Frisch \cite Scalmani2010-tw
   *
   * We implement it in terms of the EQ partition, where the weights are the
   * same for all points and there is no Gaussian smearing of the point charges:
   *
   * \f[
   *   f_{i} = 4\pi -  w_{0}\sum_{j \neq i} \frac{1}{r_{ij}},
   * \f]
   *
   * here \f$r_{ij} = |\mathbf{r}_{i} - \mathbf{r}_{j} |\f$ are the distances
   * between EQ points on the unit sphere.
   * The self-potential factors for the non-unit sphere are: \f$ f_{i}(R) =
   * Rf_{i}\f$.
   */
  Eigen::VectorXd self_potentials_;

  /** Self-field factors **for the unit sphere**.
   *
   * The formula is given by Scalmani and Frisch \cite Scalmani2010-tw It is
   * derived from Gauss' theorem and was first published by Purisima and Nilar
   * \cite Purisima1995-kc
   *
   * We implement it in terms of the EQ partition, where the weights are the
   * same for all points and there is no Gaussian smearing of the point charges:
   *
   * \f[
   *   g_{i}
   *   =
   *   - 2\pi
   *   - w_{0}\sum_{j \neq i} \frac{(\mathbf{r}_i - \mathbf{r}_{j})\cdot
   * \hat{\mathbf{n}}_{j}}{r^{3}_{ij}},
   * \f]
   *
   * here \f$r_{ij} = |\mathbf{r}_{i} - \mathbf{r}_{j} |\f$ is the distance
   * between EQ points \f$i\f$ and \f$j\f$ on the unit sphere, while
   * \f$\hat{\mathbf{n}}_{i}\f$ is the normal vector at EQ point \f$i\f$. Since
   * we are on the unit sphere at the origin, the EQ points **are** also the
   * normal vectors. The self-field factors for the non-unit sphere are: \f$
   * g_{i}(R) = g_{i}\f$.
   */
  Eigen::VectorXd self_fields_;

public:
  /** Constructor from maximum quadrature weight.
   *
   * @param[in] max_w maximum quadrature weight (inverse of the density of
   * points)
   * @param[in] r radius of the sphere.
   * @param[in] c center of the sphere.
   */
  Sphere(double max_w, double r, const Eigen::Vector3d& c, bool on_atom = true);

  /** Value of switching function at given point.
   *
   * @param[in] p evaluation point.
   * @param[in] rho interaction radius of point.
   *
   * The switching function is a smoothed Heaviside step function,
   * 0 if the point is inside this sphere, 1 if it's outside.
   *
   * Formally:
   *
   * \f[
   *   \tilde{H}(x) = \frac{1}{2}\left[ 1 + \erf\left(x\right) \right],
   * \f]
   *
   * where \f$x\f$ is the *penetration distance*:
   *
   * \f[
   *   x = \frac{s_{\alpha}(\mathbf{p})}{\rho},
   * \f]
   *
   * defined in terms of the *signed normal distance* \f$s_{\alpha}(\mathbf{p})
   * = |\mathbf{p} - \mathbf{r}_{\alpha}| - R_{\alpha}\f$ and the *interaction
   * radius* \f$\rho\f$. The latter is a function of the given evaluation point.
   */
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
  auto weight() const -> double { return w_0_ * std::pow(radius_, 2); }

  /** Weight of EQ points **on the unit sphere** */
  auto omega() const -> double { return w_0_; }

  /** EQ points. */
  auto points() const -> Eigen::Matrix3Xd { return points_; }

  auto fs() const -> Eigen::VectorXd { return radius_ * self_potentials_; }
  auto fs(size_t i) const -> double { return radius_ * self_potentials_(i); }
  auto gs() const -> Eigen::VectorXd { return self_fields_; }
  auto gs(size_t i) const -> double { return self_fields_(i); }

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
