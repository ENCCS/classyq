#include "math_utils.hpp"

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

auto surface_of_sphere(double r) -> double {
  return 4.0 * M_PI * std::pow(r, 2);
}

auto surface_of_cap(double theta, double r) -> double {
  return 2.0 * M_PI * std::pow(r, 2) * (1 - std::cos(theta));
}

auto
spherical_to_cartesian(double theta, double phi, double R,
                       const Eigen::Vector3d &c)
-> Eigen::Vector3d {
  auto sin_theta = std::sin(theta);
  auto cos_theta = std::cos(theta);
  auto sin_phi = std::sin(phi);
  auto cos_phi = std::cos(phi);

  Eigen::Vector3d r(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

  return (R * r + c);
}

auto spherical_to_cartesian(const Eigen::Matrix2Xd &sph)
-> Eigen::Matrix3Xd {
  Eigen::Array2Xd cos = sph.array().cos();
  Eigen::Array2Xd sin = sph.array().sin();

  Eigen::Matrix3Xd car = Eigen::Matrix3Xd(3, sph.cols());
  car << sin.row(0) * cos.row(1), sin.row(0) * sin.row(1), cos.row(0);

  return car;
}

auto spherical_to_cartesian(const Eigen::Matrix2Xd &sph, double R,
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

auto polar_angle_cap_of_intersection(double R, double r, double d)
-> double {
  return (std::pow(R, 2) - std::pow(r, 2) + std::pow(d, 2)) / (2.0 * d * R);
}
