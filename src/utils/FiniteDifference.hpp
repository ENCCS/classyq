#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

namespace classyq {
/**
 * Fornberg, B. Generation of Finite Difference Formulas on Arbitrarily
 * Spaced Grids. Math. Comput. 1988, 51 (184), 699–699.
 * https://doi.org/10.1090/s0025-5718-1988-0935077-0.
 */
namespace utils {

struct FivePointStencil {
  static constexpr auto NPoints = 5;

  static constexpr auto points = std::array{0, 1, -1, 2, -2};

  static constexpr auto order_1 =
      std::array{0.0, 2.0 / 3.0, -2.0 / 3.0, -1.0 / 12.0, 1.0 / 12.0};

  static constexpr auto order_2 =
      std::array{-5.0 / 2.0, 4.0 / 3.0, 4.0 / 3.0, -1.0 / 12.0, -1.0 / 12.0};

  static constexpr auto order_3 =
      std::array{0.0, -1.0, 1.0, 1.0 / 2.0, -1.0 / 12.0};

  static constexpr auto order_4 = std::array{6.0, -4.0, -4.0, 1.0, 1.0};
};

struct SevenPointStencil {
  static constexpr auto NPoints = 7;

  static constexpr auto points = std::array{0, 1, -1, 2, -2, 3, -3};

  static constexpr auto order_1 =
      std::array{0.0,        3.0 / 4.0,  -3.0 / 4.0, -3.0 / 20.0,
                 3.0 / 20.0, 1.0 / 60.0, -1.0 / 60.0};

  static constexpr auto order_2 =
      std::array{-49.0 / 18.0, 3.0 / 2.0,  3.0 / 2.0, -3.0 / 20.0,
                 -3.0 / 20.0,  1.0 / 90.0, 1.0 / 90.0};

  static constexpr auto order_3 = std::array{0.0,  -13.0 / 8.0, 13.0 / 8.0, 1.0,
                                             -1.0, -1.0 / 8.0,  1.0 / 8.0};

  static constexpr auto order_4 = std::array{
      28.0 / 3.0, -13.0 / 2.0, -13.0 / 2.0, 2.0, 2.0, -1.0 / 6.0, -1.0 / 6.0};
};

struct NinePointStencil {
  static constexpr auto NPoints = 9;

  static constexpr auto points = std::array{0, 1, -1, 2, -2, 3, -3, 4, -4};

  static constexpr auto order_1 =
      std::array{0.0,         4.0 / 5.0,    -4.0 / 5.0,   -1.0 / 5.0, 1.0 / 5.0,
                 4.0 / 105.0, -4.0 / 105.0, -1.0 / 280.0, 1.0 / 280.0};

  static constexpr auto order_2 = std::array{
      -205.0 / 72.0, 8.0 / 5.0,   8.0 / 5.0,    -1.0 / 5.0,  -1.0 / 5.0,
      8.0 / 315.0,   8.0 / 315.0, -1.0 / 560.0, -1.0 / 560.0};

  static constexpr auto order_3 = std::array{
      0.0,         -61.0 / 30.0, 61.0 / 30.0, 169.0 / 120.0, -169.0 / 120.0,
      -3.0 / 10.0, 3.0 / 10.0,   7.0 / 240.0, -7.0 / 240.0};

  static constexpr auto order_4 = std::array{
      91.0 / 8.0, -122.0 / 15.0, -122.0 / 15.0, 169.0 / 60.0, 169.0 / 60.0,
      -2.0 / 5.0, -2.0 / 5.0,    7.0 / 240.0,   7.0 / 240.0

  };
};

namespace detail {
template <typename Stencil, size_t Order,
          typename Retval = std::array<
              std::tuple<typename decltype(Stencil::points)::value_type,
                         typename decltype(Stencil::order_1)::value_type>,
              Stencil::NPoints>>
inline auto get_stencil() -> Retval {
  static_assert(Order >= 1 && Order <= 4, "Derivative order must be >=1, <=4");

  Retval retval;
  if constexpr (Order == 1) {
    std::transform(
        Stencil::points.cbegin(), Stencil::points.cend(),
        Stencil::order_1.cbegin(), retval.begin(),
        [](const auto &x, const auto &w) { return std::make_tuple(x, w); });
  } else if constexpr (Order == 2) {
    std::transform(
        Stencil::points.cbegin(), Stencil::points.cend(),
        Stencil::order_2.cbegin(), retval.begin(),
        [](const auto &x, const auto &w) { return std::make_tuple(x, w); });
  } else if constexpr (Order == 3) {
    std::transform(
        Stencil::points.cbegin(), Stencil::points.cend(),
        Stencil::order_3.cbegin(), retval.begin(),
        [](const auto &x, const auto &w) { return std::make_tuple(x, w); });
  } else {
    std::transform(
        Stencil::points.cbegin(), Stencil::points.cend(),
        Stencil::order_4.cbegin(), retval.begin(),
        [](const auto &x, const auto &w) { return std::make_tuple(x, w); });
  }

  return retval;
}
} // namespace detail

/** Finite-difference driver
 */
template <size_t Order, typename Stencil> struct FiniteDifference {
  template <typename Fn,
            typename Retval = typename std::result_of<Fn(double)>::type>
  static auto compute(Fn &&f, double h = 0.001) -> Retval {
    auto stencil = detail::get_stencil<Stencil, Order>();

    auto h_order = std::pow(h, Order);

    std::vector<Retval> tmps;
    for (const auto &[x, w] : stencil) {
      tmps.push_back((w / h_order) * f(h * x));
    }

    return std::accumulate(std::next(tmps.cbegin()), tmps.cend(), tmps[0]);
  }
};
} // namespace utils
} // namespace classyq
