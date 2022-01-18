#pragma once

#include <cmath>
#include <numeric>
#include <tuple>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "LeopardiPartition.hpp"
#include "math_utils.hpp"

/**
 */
class Sphere final {
private:
  /** Number of EQ points. */
  size_t N_{0};
  /** Weight of EQ points. */
  double w_0_;
  /** Radius. */
  double radius_;
  /** Center. */
  Eigen::Vector3d center_;
  /** EQ points, in Cartesian coordinates. */
  Eigen::Matrix3Xd points_;

public:
  /** Constructor from maximum quadrature weight.
   *
   * @param[in] max_w maximum quadrature weight.
   * @param[in] r radius of the sphere.
   * @param[in] c center of the sphere.
   */
  Sphere(double max_w, double r, Eigen::Vector3d c);

  /** Value of switching function for given point.
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
  auto switching(const Eigen::Vector3d &p, double rho) const -> double;

  /** Center of sphere. */
  auto center() const -> Eigen::Vector3d { return center_; }

  /** Radius of sphere. */
  auto radius() const -> double { return radius_; }

  /** Number of quadrature points. */
  auto size() const -> size_t { return N_; }

  /** Weight of EQ points. */
  auto w_0() const -> double { return w_0_; }

  /** EQ points. */
  auto points() const -> Eigen::Matrix3Xd { return points_; }
};

using NeighborsList = std::vector<std::vector<size_t>>;

inline auto neighbors_list(const std::vector<Sphere> &spheres)
    -> NeighborsList {
  auto ns = NeighborsList(spheres.size());

  for (auto I = 0; I < spheres.size(); ++I) {
    auto cI = spheres[I].center();
    auto rI = spheres[I].radius();
    for (auto J = 0; J < spheres.size(); ++J) {
      if (I != J) {
        auto d = (cI - spheres[J].center()).norm();
        auto r = rI + spheres[J].radius();
        if (d <= r) {
          ns[I].push_back(J);
        }
      }
    }
  }

  return ns;
}

// TODO define iterator that returns weight and point
// TODO CTOR from list of spheres
struct Quadrature final {
  size_t N_{0};
  Eigen::Matrix3Xd points_;
  Eigen::VectorXd weights_;
};

inline auto combine(const std::vector<Sphere> &spheres,
                    double threshold = std::numeric_limits<double>::epsilon())
    -> Eigen::VectorXd {

  // compute neighbors' list
  auto ns = neighbors_list(spheres);

  std::vector<double> ws;
  std::vector<double> ps;

  // loop on spheres
  for (auto I = 0; I < spheres.size(); ++I) {
    auto w_0 = spheres[I].w_0();
    auto rho = std::sqrt(w_0 / M_PI);
    // loop on sampling points on current sphere
    const auto ps_I = spheres[I].points();
    for (const auto &p_I : ps_I.colwise()) {
      auto w_I = w_0;
      // loop on interesecting spheres
      // O(N_sph^2) version
      // for (auto J = 0; J < spheres.size(); ++J) {
      //  exit early if self-intersecting
      // if (I != J) {
      // w_I *= spheres[J].switching(p_I, rho);
      //}
      // using the neighbor list, should be O(N_sph log(N_sph))
      for (const auto J : ns[I]) {
        w_I *= spheres[J].switching(p_I, rho);
      }
      if (w_I >= threshold) {
        // update counter for number of points sphere I
        ws.push_back(w_I);
        ps.push_back(p_I[0]);
        ps.push_back(p_I[1]);
        ps.push_back(p_I[2]);
      }
    }
  }

  Eigen::Matrix3Xd points =
      Eigen::Map<Eigen::Matrix3Xd>(ps.data(), 3, ps.size() / 3);

  return Eigen::Map<Eigen::VectorXd>(ws.data(), ws.size());
}
