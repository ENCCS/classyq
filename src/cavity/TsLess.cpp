#include "TsLess.hpp"

#include <algorithm>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include "Sphere.hpp"
#include "utils/utils.hpp"

namespace classyq {
TsLess::TsLess(const std::vector<Sphere>& ss, double threshold)
  : threshold_{ threshold }
  , spheres_{ ss }
{
  // compute neighbors' list
  auto neighbors = neighbors_list(ss);
  SPDLOG_TRACE("neighbors {}", neighbors);

  // find maximum number of EQ points
  auto max_EQ_points =
    (std::max_element(
       ss.cbegin(), ss.cend(), [](const Sphere& s_A, const Sphere& s_B) { return (s_A.size() < s_B.size()); }))
      ->size();
  auto max_sz = max_EQ_points * ss.size();

  // temporaries to store quadrature points, weights, normals
  std::vector<double> xs, ws, ns;
  // allocate enough space beforehand
  xs.reserve(3 * max_sz);
  ws.reserve(max_sz);
  ns.reserve(3 * max_sz);
  SPDLOG_TRACE("allocated {} of scratch for points, weights, and normal vectors",
               memory_with_units<double>(7 * max_sz));

  // loop on spheres
  for (auto A = 0; A < ss.size(); ++A) {
    auto n_on_sphere = 0;
    auto a_on_sphere = 0.0;

    auto R_A = ss[A].radius();
    auto c_A = ss[A].center();

    // EQ weight of sphere A
    auto w = ss[A].weight();
    // interaction radius of sphere A
    auto rho_A = std::sqrt(w / M_PI);

    // EQ points on sphere A a.k.a. normal vectors at quadrature points
    const auto ns_j = ss[A].points();
    // loop over EQ points on sphere A
    for (const auto& n_Aj : ns_j.colwise()) {
      // affine transformation of EQ point on sphere A a.k.a. quadrature point
      auto s_Aj = c_A + R_A * n_Aj;

      auto w_Aj = w;

      // loop on neighboring spheres
      for (const auto B : neighbors[A]) {
        w_Aj *= ss[B].switching(s_Aj, rho_A);
      }

      // accumulate quadrature points and per-sphere statistics:
      // - number of points per sphere
      // - exposed area per sphere
      if (w_Aj >= threshold) {
        n_on_sphere += 1;
        ws.push_back(w_Aj);
        a_on_sphere += w_Aj;
        // push back x, y, z of point
        xs.push_back(s_Aj(0));
        xs.push_back(s_Aj(1));
        xs.push_back(s_Aj(2));
        // push back x, y, z of normal vectors at point
        ns.push_back(n_Aj(0));
        ns.push_back(n_Aj(1));
        ns.push_back(n_Aj(2));
      }
    }
    points_per_sphere_.push_back(n_on_sphere);
    exposed_area_.push_back(a_on_sphere);
  }

  N_ = ws.size();
  weights_ = Eigen::Map<Eigen::VectorXd>(ws.data(), N_);
  points_ = Eigen::Map<Eigen::Matrix3Xd>(xs.data(), 3, N_);
  normals_ = Eigen::Map<Eigen::Matrix3Xd>(ns.data(), 3, N_);
}

auto
TsLess::str() const -> std::string
{
  // FIXME add citation
  auto str = fmt::format("TsLess cavity\n");
  // FIXME surface type
  str += fmt::format("Surface type: {}\n", "SAS");
  // FIXME
  str += fmt::format("Solvent probe radius: {:.4f}\n", 1.1);
  // FIXME conversion factor
  str += fmt::format("Average grid weight {:.4f} Ang^2\n", weights_.mean());

  auto spheres_on_atoms =
    std::count_if(spheres_.cbegin(), spheres_.cend(), [](const Sphere& s) { return s.is_atom_centered(); });
  str += fmt::format("Number of spheres: {}, of which {} centered on atoms\n", spheres_.size(), spheres_on_atoms);
  // FIXME conversion factor
  str += fmt::format("Cavity surface area: {:.4f} Ang^2\n", area());
  str += fmt::format("Number of grid points: {}\n", N_);
  // FIXME
  str += fmt::format("Number of irreducible grid points: {}\n", N_);
  str += fmt::format("============ Spheres list (in Angstrom)\n");
  // FIXME alignment and spacing
  str +=
    fmt::format("{1:s} {2:s} {3:s} ({0:s}) {4:s} {5:s} ({0:s}) {6:s} ({0:s}) {7:s} ({0:s}) {8:s} {9:s} ({0:s}^2)\n",
                "Ang",
                "Sphere",
                "on atom",
                "Radius",
                "Alpha",
                "X",
                "Y",
                "Z",
                "Grid points",
                "Surface");
  str += fmt::format("{:-^90}\n", "");
  auto idx = 0;
  // FIXME alignment and spacing
  for (const auto& s : spheres_) {
    auto c = s.center();
    // FIXME conversion factors for radius, center XYZ, exposed area
    // FIXME atom to which sphere belongs
    str += fmt::format("{0:>6d} {1:>6d} {2:>6.4f} {3:>6.2f} {4:> 6.4f} {5:> 6.4f} {6:> 6.4f} {7:>6d} {8:>6.4f}\n",
                       idx + 1,
                       0, // FIXME is not atom-centered print N/A otherwise element symbol
                       s.radius(),
                       1.0,
                       c(0),
                       c(1),
                       c(2),
                       points_per_sphere_[idx],
                       exposed_area_[idx]);
    ++idx;
  }

  return str;
}
} // namespace classyq
