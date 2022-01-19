#pragma once

#include <numeric>

#include <Eigen/Core>

#include "Sphere.hpp"

namespace classyq {
class TsLess final {
private:
  size_t N_{0};
  Eigen::VectorXd weights_;
  Eigen::Matrix3Xd points_;

public:
  TsLess(const std::vector<Sphere> &ss,
         double threshold = std::numeric_limits<double>::epsilon());

  auto size() const -> size_t { return N_; }

  auto weights() const -> Eigen::VectorXd { return weights_; }

  auto points() const -> Eigen::Matrix3Xd { return points_; }

  auto area() const -> double { return weights_.sum(); }

  /**
   * @{ Iterators
   */
  /**@}*/
};

using NeighborsList = std::vector<std::vector<size_t>>;

auto neighbors_list(const std::vector<Sphere> &spheres) -> NeighborsList;
} // namespace classyq
