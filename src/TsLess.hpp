#pragma once

#include <numeric>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Sphere.hpp"

namespace classyq {
using NeighborsList = std::vector<std::vector<size_t>>;

auto neighbors_list(const std::vector<Sphere> &spheres) -> NeighborsList;

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

template <typename T>
auto neighbors(T R0, const Eigen::Vector<T, 3> &c0, T R1,
               const Eigen::Vector<T, 3> &c1) -> bool {
  return ((c1 - c0).norm() <= R0 + R1);
}

template <typename T> auto switching(T x) -> T {
  if constexpr (std::is_same_v<T, double>) {
    return 0.5 * (1.0 + std::erf(x));
  } else {
    return 0.5 * (1.0 + erf(x));
  }
}

template <typename T>
// auto generate(const std::vector<T> &Rs,
//               const std::vector<Eigen::Vector<T, 3>> &cs, double max_w,
//               double threshold = std::numeric_limits<double>::epsilon())
auto generate(const Eigen::Vector<T, Eigen::Dynamic> &Rs,
              const Eigen::Matrix<T, 3, Eigen::Dynamic> &cs, double max_w)
    -> Eigen::Vector<T, Eigen::Dynamic> {
  // std::vector<Eigen::Vector<T, Eigen::Dynamic>> xs;
  std::vector<T> ws;

  auto N_sph = Rs.size();

  for (auto I = 0; I < N_sph; ++I) {
    T R_I = Rs(I);
    Eigen::Vector<T, 3> c_I = cs.col(I);

    // generate Leopardi partition for sphere I
    auto area = surface_of_sphere(autodiff::val(R_I));
    // number of points in Leopardi partitioning of the unit sphere.
    auto N = static_cast<size_t>(std::ceil(area / max_w));

    // Leopardi partition of the unit sphere at the origin.
    auto [w_0, points_sph] = leopardi_partition(N);

    // EQ points on the unit sphere in Cartesian coordinates
    Eigen::Matrix3Xd ps_I = spherical_to_cartesian(points_sph);

    auto w = w_0 * pow(R_I, 2);
    auto rho_Ik = sqrt(w / M_PI);

    // loop over EQ points on sphere I
    for (const auto &p : ps_I.colwise()) {
      T w_Ik = w;
      Eigen::Vector<T, 3> p_Ik = R_I * p + c_I;

      // loop over intersecting spheres
      for (auto J = 0; J < N_sph; ++J) {
        auto R_J = Rs(J);
        auto c_J = cs.col(J);

        auto do_J = (I != J) && (neighbors<T>(R_I, c_I, R_J, c_J));
        if (do_J) {

          // penetration distance of k-th point on sphere I with respect to
          // sphere J
          T x_IJk = ((p_Ik - c_J).norm() - R_J) / rho_Ik;
          w_Ik *= switching<T>(x_IJk);
        }
      }

      ws.push_back(w_Ik);
      // xs.push_back(p_Ik);
    }
  }

  Eigen::Vector<T, Eigen::Dynamic> weights =
      Eigen::Map<Eigen::Vector<T, Eigen::Dynamic>>(ws.data(), ws.size());

  return weights;
}
} // namespace classyq
