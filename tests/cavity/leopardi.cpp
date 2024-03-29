#include <doctest/doctest.h>

#include <Eigen/Core>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <highfive/H5Easy.hpp>
#include <spdlog/spdlog.h>

#include "cavity/LeopardiPartition.hpp"
#include "utils/math_utils.hpp"

#include "../utils.hpp"

using namespace classyq;

TEST_CASE("Leopardi partition with 256 points")
{
  size_t N = 1 << 8;

  auto ret = leopardi_partition(N);

  auto w_0 = std::get<0>(ret);
  REQUIRE(w_0 == doctest::Approx((4.0 * M_PI) / N));

  auto points = std::get<1>(ret);

  // load reference data from Matlab implementation
  H5Easy::File file("cavity/leopardi.h5", H5Easy::File::ReadOnly);

  // data is stored in mathematicians' convention:
  // column 0: azimuthal angles
  // column 1: polar angles
  auto ref_sph = H5Easy::load<Eigen::MatrixX2d>(file, "/leopardi/spherical/256");

  REQUIRE(allclose(points.row(0).transpose(), ref_sph.col(1)));

  REQUIRE(allclose(points.row(1).transpose(), ref_sph.col(0)));

  // column 0: x
  // column 1: y
  // column 2: z
  auto ref_cart = H5Easy::load<Eigen::MatrixX3d>(file, "/leopardi/cartesian/256");

  // perform spherical-to-Cartesian transformation
  for (auto i = 0; i < N; ++i) {
    auto p = spherical_to_cartesian(points(0, i), points(1, i));
    REQUIRE(allclose(p.transpose(), ref_cart.row(i)));
  }
}
