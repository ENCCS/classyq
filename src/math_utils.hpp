#pragma once

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

/** Surface area of a sphere.
 *
 * @param[in] r radius.
 * @return surface of the sphere \f$A = 4\pi r^{2}\f$.
 */
inline auto surface_of_sphere(double r) -> double {
  return 4.0 * M_PI * std::pow(r, 2);
}

/** Surface area of spherical cap.
 *
 * @param[in] theta polar angle between the rays from the center of the sphere
 * to the apex of the cap (the pole) and the edge of the disk forming the base
 * of the cap.
 * @param[in] r radius of the sphere.
 * @return surface of the spherical cap \f$ A = 2\pi r^{2}(1 - \cos \theta)\f$.
 *
 * @note Leopardi uses the equivalent formula:
 * \f$ A = 4 \pi r^{2} \sin^{2}\left(\frac{\theta}{2}\right)\f$.
 */
inline auto surface_of_cap(double theta, double r) -> double {
  return 2.0 * M_PI * std::pow(r, 2) * (1 - std::cos(theta));
}

/** Spherical to Cartesian transformation.
 *
 * @param[in] theta polar angle of point.
 * @param[in] phi azimuthal angle of point.
 * @param[in] R radial distance of point.
 * @param[in] c translation vector.
 */
inline auto
spherical_to_cartesian(double theta, double phi, double R = 1.0,
                       const Eigen::Vector3d &c = Eigen::Vector3d::Zero())
    -> Eigen::Vector3d {
  auto sin_theta = std::sin(theta);
  auto cos_theta = std::cos(theta);
  auto sin_phi = std::sin(phi);
  auto cos_phi = std::cos(phi);

  Eigen::Vector3d r(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

  return (R * r + c);
}

/** Spherical to Cartesian transformation on the unit sphere.
 *
 * @param[in] sph (2 x N) matrix of spherical coordinates.
 * @param[in] R radial distance of point.
 * @param[in] c translation vector.
 *
 * See here: https://stackoverflow.com/a/51569396/2528668
 */
inline auto spherical_to_cartesian(const Eigen::Matrix2Xd &sph)
    -> Eigen::Matrix3Xd {
  Eigen::Array2Xd cos = sph.array().cos();
  Eigen::Array2Xd sin = sph.array().sin();

  Eigen::Matrix3Xd car = Eigen::Matrix3Xd(3, sph.cols());
  car << sin.row(0) * cos.row(1), sin.row(0) * sin.row(1), cos.row(0);

  return car;
}

/** Spherical to Cartesian transformation.
 *
 * @param[in] sph (2 x N) matrix of spherical coordinates.
 * @param[in] R radial distance of point.
 * @param[in] c translation vector.
 *
 * See here: https://stackoverflow.com/a/51569396/2528668
 */
inline auto spherical_to_cartesian(const Eigen::Matrix2Xd &sph, double R,
                                   const Eigen::Vector3d &c)
    -> Eigen::Matrix3Xd {
  Eigen::Array2Xd cos = sph.array().cos();
  Eigen::Array2Xd sin = sph.array().sin();

  Eigen::Matrix3Xd car = Eigen::Matrix3Xd(3, sph.cols());
  car << sin.row(0) * cos.row(1), sin.row(0) * sin.row(1), cos.row(0);

  Eigen::Transform<double, 3, Eigen::Affine> T =
      Eigen::Translation3d(c) * Eigen::Scaling(R);

  return (T * car);
}

/** Colatitude of the spherical cap resulting from the intersection of two
 * spheres.
 *
 * @param[in] R radius of sphere 1.
 * @param[in] r radius of sphere 2.
 * @param[in] d distance between the two spheres.
 *
 * This function returns the polar angle of the cap on **sphere 1**. It's enough
 * to exchange the radii in input to obtain the polar angle of the cap on
 * sphere 2.
 * For the derivation see:
 * https://en.wikipedia.org/wiki/Spherical_cap#Areas_of_intersecting_spheres
 */
inline auto polar_angle_cap_of_intersection(double R, double r, double d)
    -> double {
  return (std::pow(R, 2) - std::pow(r, 2) + std::pow(d, 2)) / (2.0 * d * R);
}
