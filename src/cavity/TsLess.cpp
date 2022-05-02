#include "TsLess.hpp"

#include <tuple>
#include <vector>

#include <Eigen/Core>
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/dual/eigen.hpp"
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include "Sphere.hpp"

namespace classyq {
TsLess::TsLess(const std::vector<Sphere> &spheres, double threshold) {
  // compute neighbors' list
  auto neighbors = neighbors_list(spheres);

  std::vector<double> ws;
  std::vector<double> ps;
  std::vector<double> fs;
  std::vector<double> gs;
  std::vector<double> ns;

  // loop on spheres
  for (auto I = 0; I < spheres.size(); ++I) {
    // affine transformation from unit sphere
    Eigen::Transform<double, 3, Eigen::Affine> T =
        Eigen::Translation3d(spheres[I].center()) *
        Eigen::Scaling(spheres[I].radius());

    // number of points on sphere I
    auto n_on_sphere = 0;
    // get EQ weight (same for all point on sphere I)
    auto w = spheres[I].weight();
    weights_0_.push_back(spheres[I].weight0());
    radii_.push_back(spheres[I].radius());

    auto rho = std::sqrt(w / M_PI);
    // normal vectors
    const auto ns_I = spheres[I].points();
    // loop on sampling points on current sphere
    auto k = 0;
    for (const auto &n_I : ns_I.colwise()) {
      // affine transformation
      auto p_I = T * n_I;
      autodiff::dual w_I = w;
      for (const auto J : neighbors[I]) {
        w_I *= spheres[J].switching(p_I, rho);
      }
      //  accumulate points and per-sphere statistics (e.g. exposed area, number
      //  of points)
      if (w_I >= threshold) {
        ws.push_back(autodiff::val(w_I));
        // push back x, y, z of point
        ps.push_back(p_I(0));
        ps.push_back(p_I(1));
        ps.push_back(p_I(2));
        fs.push_back(spheres[I].fs(k));
        gs.push_back(spheres[I].gs(k));
        // push back x, y, z of normal vectors at point
        ns.push_back(n_I(0));
        ns.push_back(n_I(1));
        ns.push_back(n_I(2));
        n_on_sphere += 1;
      }
      k += 1;
    }
    points_per_sphere_.push_back(n_on_sphere);
  }

  N_ = ws.size();
  weights_ = Eigen::Map<Eigen::VectorXd>(ws.data(), ws.size());
  points_ = Eigen::Map<Eigen::Matrix3Xd>(ps.data(), 3, ps.size() / 3);
  self_potentials_ = Eigen::Map<Eigen::VectorXd>(fs.data(), fs.size());
  self_fields_ = Eigen::Map<Eigen::VectorXd>(gs.data(), gs.size());
  normals_ = Eigen::Map<Eigen::Matrix3Xd>(ns.data(), 3, ns.size() / 3);

  SPDLOG_INFO("npoints_on_sphere_ {}", points_per_sphere_);
}
} // namespace classyq
