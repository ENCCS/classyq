#include <doctest/doctest.h>

#include <Eigen/Core>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <highfive/H5Easy.hpp>

#include "cavity/Sphere.hpp"
#include "cavity/TsLess.hpp"
#include "utils/math_utils.hpp"

#include "../utils.hpp"

using namespace classyq;
using namespace doctest;

TEST_CASE("TsLess partition of two interlocking spheres") {
  auto max_w = 0.1;
  auto d = 2.5;

  Eigen::Vector3d c1(0.0, 0.0, 0.0), c2(0.0, 0.0, d);
  auto r1 = 3.0;
  Sphere s1(max_w, r1, c1);

  auto r2 = 1.5;
  Sphere s2(max_w, r2, c2);

  auto a1 = surface_of_sphere(r1);
  auto cos1 = polar_angle_cap_of_intersection(r1, r2, d);
  auto cap1 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1);

  auto a2 = surface_of_sphere(r2);
  auto cos2 = polar_angle_cap_of_intersection(r2, r1, d);
  auto cap2 = 2.0 * M_PI * std::pow(r2, 2) * (1.0 - cos2);

  auto ref_area = (a1 + a2) - (cap1 + cap2);
  fmt::print("ref_area = {}\n", ref_area);

  auto tsless = TsLess({s1, s2}, 1.0e-9);

  auto area = tsless.area();
  fmt::print("area = {}\n", area);

  REQUIRE(area == Approx(ref_area).epsilon(max_w));
}

TEST_CASE("TsLess partition of three interlocking spheres") {
  // three spheres in an isosceles triangle configuration
  auto max_w = 0.1;
  // distance between centers of the interlocking spheres
  auto d = 3.0;

  Eigen::Vector3d c1(0.0, 0.0, 0.0);
  auto r1 = 2.0;
  Sphere s1(max_w, r1, c1);

  Eigen::Vector3d c2(0.0, -d * std::cos(M_PI / 4), -d * std::sin(M_PI / 4));
  auto r2 = 1.5;
  Sphere s2(max_w, r2, c2);

  Eigen::Vector3d c3(0.0, d * std::cos(M_PI / 4), -d * std::sin(M_PI / 4));
  auto r3 = 1.8;
  Sphere s3(max_w, r3, c3);

  auto a1 = surface_of_sphere(r1);
  // s1 interlocks with s2 and s3
  // cosine of polar angle for interlock of s1 and s2
  auto cos1_2 = polar_angle_cap_of_intersection(r1, r2, d);
  // surface area of cap with polar angle cos1_2
  auto cap1_2 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1_2);
  // cosine of polar angle for interlock of s1 and s3
  auto cos1_3 = polar_angle_cap_of_intersection(r1, r3, d);
  // surface area of cap with polar angle cos1_3
  auto cap1_3 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1_3);

  auto a2 = surface_of_sphere(r2);
  // s2 interlocks with s1 only
  // cosine of polar angle for interlock of s2 and s1
  auto cos2_1 = polar_angle_cap_of_intersection(r2, r1, d);
  // surface area of cap with polar angle cos2_1
  auto cap2_1 = 2.0 * M_PI * std::pow(r2, 2) * (1.0 - cos2_1);

  auto a3 = surface_of_sphere(r3);
  // s3 interlocks with s1 only
  // cosine of polar angle for interlock of s3 and s1
  auto cos3_1 = polar_angle_cap_of_intersection(r3, r1, d);
  // surface area of cap with polar angle cos3_1
  auto cap3_1 = 2.0 * M_PI * std::pow(r3, 2) * (1.0 - cos3_1);

  auto ref_area = (a1 + a2 + a3) - (cap1_2 + cap1_3 + cap2_1 + cap3_1);
  fmt::print("ref_area = {}\n", ref_area);

  auto tsless = TsLess({s1, s2, s3}, 1.0e-9);

  auto area = tsless.area();
  fmt::print("area = {}\n", area);

  REQUIRE(area == Approx(ref_area).epsilon(max_w));
}

TEST_CASE("TsLess partition of four interlocking spheres") {
  // four spheres: as above (isosceles triangle configuration) plus one sphere
  // along the z axis
  auto max_w = 0.1;
  // distance between centers of the interlocking spheres
  auto d = 3.0;

  Eigen::Vector3d c1(0.0, 0.0, 0.0);
  auto r1 = 2.0;
  Sphere s1(max_w, r1, c1);

  Eigen::Vector3d c2(0.0, -d * std::cos(M_PI / 4), -d * std::sin(M_PI / 4));
  auto r2 = 1.5;
  Sphere s2(max_w, r2, c2);

  Eigen::Vector3d c3(0.0, d * std::cos(M_PI / 4), -d * std::sin(M_PI / 4));
  auto r3 = 1.8;
  Sphere s3(max_w, r3, c3);

  Eigen::Vector3d c4(0.0, 0.0, 2 * d);
  auto r4 = 6.0;
  Sphere s4(max_w, r4, c4);

  auto a1 = surface_of_sphere(r1);
  // s1 interlocks with s2, s3, and s4
  // cosine of polar angle for interlock of s1 and s2
  auto cos1_2 = polar_angle_cap_of_intersection(r1, r2, d);
  // surface area of cap with polar angle cos1_2
  auto cap1_2 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1_2);
  // cosine of polar angle for interlock of s1 and s3
  auto cos1_3 = polar_angle_cap_of_intersection(r1, r3, d);
  // surface area of cap with polar angle cos1_3
  auto cap1_3 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1_3);
  // cosine of polar angle for interlock of s1 and s4
  auto cos1_4 = polar_angle_cap_of_intersection(r1, r4, 2 * d);
  // surface area of cap with polar angle cos1_3
  auto cap1_4 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1_4);

  auto a2 = surface_of_sphere(r2);
  // s2 interlocks with s1 only
  // cosine of polar angle for interlock of s2 and s1
  auto cos2_1 = polar_angle_cap_of_intersection(r2, r1, d);
  // surface area of cap with polar angle cos2_1
  auto cap2_1 = 2.0 * M_PI * std::pow(r2, 2) * (1.0 - cos2_1);

  auto a3 = surface_of_sphere(r3);
  // s3 interlocks with s1 only
  // cosine of polar angle for interlock of s3 and s1
  auto cos3_1 = polar_angle_cap_of_intersection(r3, r1, d);
  // surface area of cap with polar angle cos3_1
  auto cap3_1 = 2.0 * M_PI * std::pow(r3, 2) * (1.0 - cos3_1);

  auto a4 = surface_of_sphere(r4);
  // s4 interlocks with s1 only
  // cosine of polar angle for interlock of s3 and s1
  auto cos4_1 = polar_angle_cap_of_intersection(r4, r1, 2 * d);
  // surface area of cap with polar angle cos4_1
  auto cap4_1 = 2.0 * M_PI * std::pow(r4, 2) * (1.0 - cos4_1);

  auto ref_area = (a1 + a2 + a3 + a4) -
                  (cap1_2 + cap1_3 + cap1_4 + cap2_1 + cap3_1 + cap4_1);
  fmt::print("ref_area = {}\n", ref_area);

  auto tsless = TsLess({s1, s2, s3, s4}, 1.0e-13);

  auto area = tsless.area();
  fmt::print("area = {}\n", area);

  REQUIRE(area == Approx(ref_area).epsilon(max_w));
}
