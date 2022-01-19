#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <vector>

#include <Eigen/Core>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include "Sphere.hpp"
#include "TsLess.hpp"

using namespace classyq;

int main() {
  spdlog::set_pattern("[%Y-%m-%d %T][%^%l%$][TID: %t, PID: %P][%!@%s:%4#] %v");

  auto max_w = 0.1;
  auto d = 2.5;

  Eigen::Vector3d c1, c2;
  c1 << 0.0, 0.0, 0.0;
  Sphere s1(max_w, 3.0, c1);

  c2 << 0.0, 0.0, d;
  Sphere s2(max_w, 1.5, c2);

  auto cos_theta = [](double R, double r, double d) {
    return (std::pow(R, 2) - std::pow(r, 2) + std::pow(d, 2)) / (2.0 * d * R);
  };
  auto cos1 = cos_theta(s1.radius(), s2.radius(), d);
  auto cap1 = 2.0 * M_PI * std::pow(s1.radius(), 2) * (1.0 - cos1);
  auto a1 = 4.0 * M_PI * std::pow(s1.radius(), 2);
  SPDLOG_INFO("surface of s1 {}", a1);

  auto cos2 = cos_theta(s2.radius(), s1.radius(), d);
  auto cap2 = 2.0 * M_PI * std::pow(s2.radius(), 2) * (1.0 - cos2);
  auto a2 = 4.0 * M_PI * std::pow(s2.radius(), 2);
  SPDLOG_INFO("surface of s2 {}", a2);

  SPDLOG_INFO("surface of union {}", a1 - cap1 + a2 - cap2);

  auto tsless = TsLess({s1, s2});

  SPDLOG_INFO("ws.size() {}", tsless.size());
  SPDLOG_INFO("ws.sum() {}", tsless.area());

  return EXIT_SUCCESS;
}
