#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <vector>

#include <Eigen/Core>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

constexpr auto SURFACE_OF_SPHERE = 4.0 * M_PI;

constexpr auto spherical_radius_of_cap(double area) -> double {
  return 2.0 * std::asin(std::sqrt(area / M_PI) / 2.0);
}

constexpr auto surface_of_cap(double spherical_radius) -> double {
  return 4.0 * M_PI * std::pow(std::sin(spherical_radius / 2.0), 2);
}

// area of spherical collar between a_top (top spherical radius) and a_bot
// (bottom spherical radius)
constexpr auto surface_of_collar(double r_top, double r_bot) -> double {
  return surface_of_cap(r_bot) - surface_of_cap(r_top);
}

// number of regions in each collar
// returns a vector of size n_collars + 2. The first and last elements are 1
// becuase north and south pole only have one region.
auto ideal_number_of_regions(double surface_ideal_region, double polar_colat,
                             size_t n_collars) -> std::vector<size_t> {
  auto retval = std::vector<size_t>(n_collars + 2, 1);
  auto a_fitting = (M_PI - 2.0 * polar_colat) / n_collars;

  for (auto i = 0; i < n_collars; ++i) {
    auto ideal_collar_area = surface_of_collar(
        polar_colat + i * a_fitting, polar_colat + (i + 1) * a_fitting);
    retval[i + 1] = static_cast<size_t>(
        std::round(ideal_collar_area / surface_ideal_region));
  }

  return retval;
}

// list of colatitudes of caps
// returns a vector of size n_collars + 2.
// - The first element (north pole) is the polar colatitude
// - Second-to-last element is M_PI minus the polar colatitude
// - Last element in M_PI
auto cap_colatitudes(double area_ideal_region, double polar_colat,
                     const std::vector<size_t> &n_regions)
    -> std::vector<double> {
  auto retval = std::vector<double>(n_regions.size(), 0.0);

  retval.front() = polar_colat;
  retval.back() = M_PI;

  auto n_collars = n_regions.size() - 2;

  auto subtotal_n_regions = 1;

  for (auto i = 1; i < n_collars + 1; ++i) {
    subtotal_n_regions += n_regions[i];
    retval[i] = spherical_radius_of_cap(subtotal_n_regions * area_ideal_region);
  }

  return retval;
}

auto leopardi_partition(size_t N) -> Eigen::MatrixX2d {

  // this will be our weight
  auto area_ideal_region = SURFACE_OF_SPHERE / N;

  // colatitude of the north pole spherical cap (spherical radius of spherical
  // cap of given area)
  // this is theta_c in the paper
  auto theta_c = 2.0 * std::asin(std::sqrt(1.0 / N));
  // the colatitude of the North pole is theta_c, the colatitude of the South
  // pole is M_PI - theta_c

  // this is delta_I in the paper
  auto ideal_collar_angle = std::sqrt(area_ideal_region);

  // IDEAL number of collars between polar caps
  // this is n_I in the paper
  auto n_collars = static_cast<size_t>(
      std::max(1.0, std::round((M_PI - 2 * theta_c) / ideal_collar_angle)));

  SPDLOG_INFO("n_collars {}", n_collars);

  // compute fitting collar angle
  // this is delta_F in the paper
  auto a_fitting = (M_PI - 2.0 * theta_c) / n_collars;

  // this is the y array in the paper
  auto n_regions_per_collar = std::vector<size_t>(n_collars + 2, 1);

  // n_regions_per_collars[0] = 1 (only one region in the North pole)
  // n_regions_per_collars[n_collars+1] = 1 (only one region in the South pole)
  // I use a different formula from the one in the paper: we specialize on the
  // 2-sphere so we can save some FLOPs.
  // we first compute the ideal number of regions (non-integer), add the
  // accumulated difference between ideal and actual, and round to integer. this
  // guarantees that we get the requested number of regions
  auto diff = 0.0;
  auto n_regions = 2;
  for (auto i = 1; i <= n_collars; ++i) {
    auto ideal = N * std::sin(theta_c + (i - 0.5) * a_fitting) *
                 std::sin(a_fitting / 2.0);
    n_regions_per_collar[i] = std::round(ideal + diff);
    diff += ideal - n_regions_per_collar[i];
    n_regions += n_regions_per_collar[i];
  }
  SPDLOG_INFO("number of regions per collar {}", n_regions_per_collar);
  // TODO check that we generate exactly as many regions as was requested
  SPDLOG_INFO("n_regions {}", n_regions);

  // this is the theta array in the paper
  auto cap_colats = std::vector<double>(n_collars + 2, 0.0);
  cap_colats.front() = theta_c;
  // FIXME shouldn't this be M_PI - theta_c
  cap_colats.back() = M_PI;

  auto acc = static_cast<double>(n_regions_per_collar[0]);
  for (auto i = 1; i <= n_collars; ++i) {
    acc += n_regions_per_collar[i];
    cap_colats[i] = 2.0 * std::asin(std::sqrt(acc / N));
  }
  SPDLOG_INFO("cap colatitudes {}", cap_colats);

  // col 0: values of theta (aka colatitudes aka polar angles)
  // col 1: values of phi (aka azimuthal angles)
  Eigen::MatrixX2d partition_sph = Eigen::MatrixX2d::Zero(N, 2);
  // polar angle for the south pole
  partition_sph(N - 1, 0) = M_PI;

  // loop over collars, excluding north and south poles
  auto start = n_regions_per_collar[0];
  auto offset = 0.0;
  for (auto i = 1; i <= n_collars; ++i) {
    // a_top is the colatitude of the top of the current collar.
    auto a_top = cap_colats[i - 1];
    // a_bot is the colatitude of the bottom of the current collar.
    auto a_bot = cap_colats[i];
    auto n_regions = n_regions_per_collar[i];

    auto a_point = (a_top + a_bot) / 2.0;

    // all regions in this collar have the same polar angle
    partition_sph(Eigen::seqN(start, n_regions), 0).setConstant(a_point);

    // each region is now a circle, which we can partition trivially
    // loop over regions
    auto inc = 2.0 * M_PI / n_regions;
    for (auto j = 1; j <= n_regions; ++j) {
      auto x = j * inc - M_PI / n_regions;
      // azimuthal angles
      partition_sph(j + (start - 1), 1) =
          std::fmod(x + 2 * M_PI * offset, 2 * M_PI);
    }

    auto n_regions_next = n_regions_per_collar[i + 1];
    offset += 0.5 * (1.0 / n_regions_next - 1.0 / n_regions) +
              std::gcd(n_regions, n_regions_next) /
                  (2.0 * n_regions * n_regions_next);
    offset -= std::floor(offset);

    start += n_regions;
  }

  return partition_sph;
}

int main() {
  spdlog::set_pattern("[%Y-%m-%d %T][%^%l%$][TID: %t, PID: %P][%!@%s:%4#] %v");

  // number of points in partition
  constexpr auto N = (1 << 5);

  auto points_sph = leopardi_partition(N);

  SPDLOG_INFO("points_sph\n{}", points_sph);

  return EXIT_SUCCESS;
}
