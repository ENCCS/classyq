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

  std::vector<double> ws;
  std::vector<double> ps;
  std::vector<double> fs;

  std::vector<double> radii;

  // loop on spheres
  for (auto I = 0; I < spheres.size(); ++I) {
    // affine transformation from unit sphere
    Eigen::Transform<double, 3, Eigen::Affine> T =
        Eigen::Translation3d(spheres[I].center()) *
        Eigen::Scaling(spheres[I].radius());

    // number of points on sphere I
    auto npoints_on_sphere = 0;
    // get EQ weight (same for all point on sphere I)
    auto w = spheres[I].weight();
    auto rho = std::sqrt(w / M_PI);
    // loop on sampling points on current sphere
    const auto ps_I = spheres[I].points();
    auto k = 0;
    for (const auto &p_I : ps_I.colwise()) {
      auto w_I = w;
      for (const auto J : ns[I]) {
        w_I *= spheres[J].switching(p_I, rho);
      }
      // accumulate points and per-sphere statistics (e.g. exposed area, number
      // of points)
      if (w_I >= threshold) {
        ws.push_back(w_I);
        // push back x, y, z of point
        ps.push_back(p_I(0));
        ps.push_back(p_I(1));
        ps.push_back(p_I(2));
        fs.push_back(spheres[I].fs(k));
        auto n_I = (p_I - spheres[I].center()) / spheres[I].radius();
        // TODO save radius of point
        radii[k] = spheres[I].radius();
        // TODO accumulate self-field factors
        npoints_on_sphere += 1;
      }
      k += 1;
    }
    SPDLOG_INFO("npoints_on_sphere {}", npoints_on_sphere);
  }

  N_ = ws.size();
  SPDLOG_INFO("N_ {}", N_);
  weights_ = Eigen::Map<Eigen::VectorXd>(ws.data(), ws.size());
  points_ = Eigen::Map<Eigen::Matrix3Xd>(ps.data(), 3, ps.size() / 3);
  self_potentials_ = Eigen::Map<Eigen::VectorXd>(fs.data(), fs.size());
  // self_fields_ = Eigen::Map<Eigen::VectorXd>(gs.data(), gs.size());
}
} // namespace classyq
