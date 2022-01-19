#include <catch2/catch.hpp>

#include <Eigen/Core>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <highfive/H5Easy.hpp>
#include <spdlog/spdlog.h>

#include "Sphere.hpp"
#include "TsLess.hpp"
#include "math_utils.hpp"

#include "../utils.hpp"

using namespace classyq;

TEST_CASE("TsLess partition of two interlocking spheres",
          "[cavity][tsless][two_spheres]") {
  auto max_w = 0.1;
  auto d = 2.5;

  Eigen::Vector3d c1(0.0, 0.0, 0.0), c2(0.0, 0.0, d);
  auto r1 = 3.0;
  Sphere s1(max_w, r1, c1);

  auto r2 = 1.5;
  Sphere s2(max_w, r2, c2);

  auto tsless = TsLess({s1, s2});

  auto area = tsless.area();

  auto a1 = surface_of_sphere(r1);
  auto cos1 = polar_angle_cap_of_intersection(r1, r2, d);
  auto cap1 = 2.0 * M_PI * std::pow(r1, 2) * (1.0 - cos1);

  auto a2 = surface_of_sphere(r2);
  auto cos2 = polar_angle_cap_of_intersection(r2, r1, d);
  auto cap2 = 2.0 * M_PI * std::pow(r2, 2) * (1.0 - cos2);

  auto ref_area = (a1 + a2) - (cap1 + cap2);

  REQUIRE(area == Approx(ref_area).epsilon(max_w));
}
