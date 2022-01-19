#include "TsLess.hpp"

#include <vector>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "Sphere.hpp"

namespace classyq {
auto neighbors_list(const std::vector<Sphere> &spheres) -> NeighborsList {
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

TsLess::TsLess(const std::vector<Sphere> &spheres, double threshold) {
  // compute neighbors' list
  auto ns = neighbors_list(spheres);

  SPDLOG_TRACE("neighbor list = {}", ns);

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

  N_ = ws.size();
  weights_ = Eigen::Map<Eigen::VectorXd>(ws.data(), ws.size());
  points_ = Eigen::Map<Eigen::Matrix3Xd>(ps.data(), 3, ps.size() / 3);
}
} // namespace classyq
