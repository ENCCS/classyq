#pragma once

#include <Eigen/Core>

namespace classyq {
/** Surface area of a sphere.
 *
 * @param[in] r radius.
 * @return surface of the sphere \f$A = 4\pi r^{2}\f$.
 */
auto
surface_of_sphere(double r) -> double;

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
auto
surface_of_cap(double theta, double r) -> double;

/** Spherical to Cartesian transformation.
 *
 * @param[in] theta polar angle of point.
 * @param[in] phi azimuthal angle of point.
 * @param[in] R radial distance of point.
 * @param[in] c translation vector.
 */
auto
spherical_to_cartesian(double theta, double phi, double R = 1.0, const Eigen::Vector3d& c = Eigen::Vector3d::Zero())
  -> Eigen::Vector3d;

/** Spherical to Cartesian transformation on the unit sphere.
 *
 * @param[in] sph (2 x N) matrix of spherical coordinates.
 *
 * See here: https://stackoverflow.com/a/51569396/2528668
 */
auto
spherical_to_cartesian(const Eigen::Matrix2Xd& sph) -> Eigen::Matrix3Xd;

/** Spherical to Cartesian transformation.
 *
 * @param[in] sph (2 x N) matrix of spherical coordinates.
 * @param[in] R radial distance of point.
 * @param[in] c translation vector.
 *
 * See here: https://stackoverflow.com/a/51569396/2528668
 */
auto
spherical_to_cartesian(const Eigen::Matrix2Xd& sph, double R, const Eigen::Vector3d& c) -> Eigen::Matrix3Xd;

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
auto
polar_angle_cap_of_intersection(double R, double r, double d) -> double;
} // namespace classyq
