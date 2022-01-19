#pragma once

#include <cmath>

#include <Eigen/Core>

#include "LeopardiPartition.hpp"
#include "math_utils.hpp"

namespace classyq {
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
} // namespace classyq
