#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include "cavity/LeopardiPartition.hpp"
#include "cavity/Sphere.hpp"
#include "cavity/TsLess.hpp"
#include "solver/Solver.hpp"
#include "solver/SolverImpl.hpp"
#include "utils/FiniteDifference.hpp"

using namespace classyq;

Eigen::VectorXd
computeMEP(const Eigen::Matrix3Xd& grid, double charge, const Eigen::Vector3d& origin)
{
  return (charge / (grid.colwise() - origin).colwise().norm().array()).colwise().sum();
}

using Stencil = std::map<int, std::vector<double>>;

// auto first = Stencil{
//     {3, std::vector{-1.0 / 2.0, 0.0, 1.0 / 2.0}},
//     {5, std::vector{1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0}},
//     {9, std::vector{1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0.0,
//                     4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0}}};

auto first = Stencil{ { 5, std::vector{ 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 } } };

auto stencils = std::map<int, Stencil>{ { 1, first } };

constexpr auto second_derivative =
  std::array<double, 9>{ -1.0 / 560.0, 8.0 / 315.0, -1.0 / 5.0,  8.0 / 5.0,   -205.0 / 72.0,
                         8.0 / 5.0,    -1.0 / 5.0,  8.0 / 315.0, -1.0 / 560.0 };

constexpr auto third_derivative = std::array<double, 9>{ -7.0 / 240.0, 3.0 / 10.0,    -169.0 / 120.0, 61.0 / 30.0, 0.0,
                                                         -61.0 / 30.0, 169.0 / 120.0, -3.0 / 10.0,    7.0 / 240.0 };

// for 2-atom molecule
auto
finite_difference(double h = 0.001)
{
  auto density = 0.3;

  Eigen::Vector3d c1(0.0, -3.0 * std::cos(M_PI / 4), -3.0 * std::sin(M_PI / 4));
  Sphere s1(density, 1.5, c1);

  // build derivative wrt center of s1
  for (const auto& [order, stencil] : stencils) {
    SPDLOG_INFO("Computing derivative of order {}", order);
    for (const auto& [npoints, coefs] : stencil) {
      SPDLOG_INFO("Using {}-point (centered) stencil", npoints);

      std::vector<Eigen::VectorXd> tmps;

      int i = -(coefs.size() / 2);

      for (const auto c : coefs) {
        Eigen::Vector3d c0(0.0, 0.0, h * i);
        Sphere s0(density, 2.0, c0);
        Eigen::VectorXd ws = TsLess({ s0, s1 }, 1.0e-9).weights();
        SPDLOG_INFO("n_points {}", ws.size());

        tmps.push_back((c / h) * ws);

        i += 1;
      }

      Eigen::VectorXd dwds1_z = std::accumulate(std::next(tmps.cbegin()), tmps.cend(), tmps[0]);

      SPDLOG_INFO("dwds1_z\n{}", dwds1_z);
    }
  }

  // derivative wrt center of s2
}

int
main()
{
  spdlog::set_pattern("[%Y-%m-%d %T][%^%l%$][TID: %t, PID: %P][%!@%s:%4#] %v");

  auto max_w = 0.1;
  auto permittivity = 4.23;

  Eigen::Vector3d c1{0.0, 0.0, 0.0};
  auto r1 = 3.0;
  Sphere s1(max_w, r1, c1);

  auto d = 2.5;

  Eigen::Vector3d c2{0.0, 0.0, d};
  auto r2 = 1.5;
  Sphere s2(max_w, r2, c2);

  auto tsless = TsLess({ s1, s2 }, 1.0e-9);
  SPDLOG_INFO("tsless {}", tsless);

  // compute S, we enforce Hermiticity
  auto S = form_S(tsless);

  // CPCM
  // generate RHS vector: potential of a point charge at c1
  auto charge = 8.0;
  auto totalASC = -charge * (permittivity - 1) / (permittivity);
  auto V = computeMEP(tsless.points(), charge, Eigen::Vector3d{ 0.0, 0.0, 0.8 });
  // solve for the charges
  Eigen::VectorXd q = -(permittivity - 1.0) / (permittivity)*S.ldlt().solve(V);
  // Gauss estimate
  auto gauss = q.sum();
  SPDLOG_INFO("<< CPCM >>\ntotalASC = {}; gauss = {}; Delta = {}", totalASC, gauss, totalASC - gauss);

  // isotropic IEFPCM
  auto R_epsilon = form_R_epsilon(tsless, permittivity);
  auto R_infinity = form_R_infinity(tsless);
  auto T = R_epsilon * S;

  // this is not using the polarization weights!
  q = -T.lu().solve(R_infinity * V);
  gauss = q.sum();
  SPDLOG_INFO("<< IEFPCM >>\ntotalASC = {}; gauss = {}; Delta = {}", totalASC, gauss, totalASC - gauss);

  //  finite_difference();
  //
  //  auto density = 0.3;
  //  auto gen = [density, c1](double z) {
  //    Sphere s1(density, 1.5, c1);
  //    Eigen::Vector3d c0(0.0, 0.0, z);
  //    Sphere s0(density, 2.0, c0);
  //    return TsLess({ s0, s1 }, 1.0e-9).weights();
  //  };
  //
  //  auto dwds0_z = utils::FiniteDifference<1, utils::FivePointStencil>::compute(gen);
  //
  //  SPDLOG_INFO("dwds0_z\n{}", dwds0_z);

  return EXIT_SUCCESS;
}
