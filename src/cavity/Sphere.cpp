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

#include "Sphere.hpp"

#include <cmath>
#include <tuple>

#include <Eigen/Core>
#include <fmt/core.h>
#include <fmt/ostream.h>

#include "LeopardiPartition.hpp"
#include "utils/math_utils.hpp"

namespace classyq {
Sphere::Sphere(double max_w, double r, const Eigen::Vector3d& c, bool on_atom)
  : atom_centered_{ on_atom }
  , radius_{ r }
  , center_{ c }
{
  auto area = surface_of_sphere(radius_);
  // number of points in Leopardi partitioning of the unit sphere.
  N_ = static_cast<size_t>(std::ceil(area / max_w));

  // Leopardi partition of the unit sphere at the origin.
  Eigen::Matrix2Xd points_sph;
  std::tie(omega_, points_sph) = leopardi_partition(N_);

  // compute normal vectors at the EQ points by transforming EQ points from
  // spherical to Cartesian coordinates on the unit sphere
  points_ = spherical_to_cartesian(points_sph);
}

auto
Sphere::switching(const Eigen::Vector3d& s, double rho) const -> double
{
  // penetration distance
  auto d = ((s - center_).norm() - radius_) / rho;

  return 0.5 * (1.0 + std::erf(d));
}

auto
Sphere::str() const -> std::string
{
  auto str = fmt::format("Center: {}\n", center_);
  // FIXME conversion factor
  str += fmt::format("Radius: {:.4f} (Ang)\n", radius_);
  // FIXME conversion factor
  str += fmt::format("Surface area: {:.4f} (Ang^2)\n", surface_of_sphere(radius_));
  str += fmt::format("Number of EQ points: {}\n", N_);
  str += fmt::format("EQ weight: {:.4f}\n", weight());
  str += fmt::format("Atom centered? {}", atom_centered_);

  return str;
}

auto
neighbors_list(const std::vector<Sphere>& ss) -> NeighborsList
{
  auto ns = NeighborsList(ss.size());

  for (auto A = 0; A < ss.size(); ++A) {
    auto c_A = ss[A].center();
    auto R_A = ss[A].radius();
    for (auto B = 0; B < ss.size(); ++B) {
      auto c_B = ss[B].center();
      auto R_B = ss[B].radius();
      if (A != B) {
        auto d = (c_A - c_B).norm();
        auto r = R_A + R_B;
        if (d <= r) {
          ns[A].push_back(B);
        }
      }
    }
  }

  return ns;
}
} // namespace classyq
