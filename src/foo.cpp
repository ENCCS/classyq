#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include "Sphere.hpp"
#include "TsLess.hpp"

using namespace classyq;

Eigen::VectorXd computeMEP(const Eigen::Matrix3Xd &grid, double charge,
                           const Eigen::Vector3d &origin) {
  return (charge / (grid.colwise() - origin).colwise().norm().array())
      .colwise()
      .sum();
}

int main() {
  spdlog::set_pattern("[%Y-%m-%d %T][%^%l%$][TID: %t, PID: %P][%!@%s:%4#] %v");

  auto max_w = 0.1;
  auto permittivity = 4.23;

  Eigen::Vector3d c1(0.0, 0.0, 0.0);
  auto r1 = 1.0;
  Sphere s1(max_w, r1, c1);

  auto d = 2.5;

  Eigen::Vector3d c2(0.0, 0.0, d);
  auto r2 = 1.5;
  Sphere s2(max_w, r2, c2);

  auto tsless = TsLess({s1, s2});

  const auto sz = tsless.size();

  // compute S, we enforce Hermiticity
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(sz, sz);
  S.diagonal() = tsless.fs().array() / tsless.weights().array();
  for (auto i = 0; i < sz; ++i) {
    const auto p_i = tsless.points(i);
    for (auto j = i + 1; j < sz; ++j) {
      S(i, j) = S(j, i) = 1.0 / (p_i - tsless.points(j)).norm();
    }
  }

  // Eigen::MatrixXd D = Eigen::MatrixXd::Zero(sz, sz);
  // D.diagonal() = tsless_1.gs().array() / tsless_1.weights().array();

  // CPCM
  // generate RHS vector: potential of a point charge at c1
  auto charge = 8.0;
  auto totalASC = -charge * (permittivity - 1) / (permittivity);
  auto V = computeMEP(tsless.points(), charge, Eigen::Vector3d(0.0, 0.0, 1.25));
  // solve for the charges
  Eigen::VectorXd q = -(permittivity - 1.0) / (permittivity)*S.ldlt().solve(V);
  // Gauss estimate
  auto gauss = q.sum();
  SPDLOG_INFO("totalASC = {}; gauss = {}; Delta = {}", totalASC, gauss,
              totalASC - gauss);

  return EXIT_SUCCESS;
}
