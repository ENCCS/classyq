#include <cmath>
#include <cstdlib>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

using namespace autodiff;

template <typename T>
auto generate(const Eigen::Vector3d &C, const Eigen::Vector<T, 3> &x) -> T {
  return 1.0 / (C - x).norm();
}

template <typename T>
auto generate2(const Eigen::Vector3d &C,
               const Eigen::Matrix<T, 3, Eigen::Dynamic> &xs) -> T {
  return 1.0 / (C - xs.rowwise().sum()).norm();
}

template <typename T>
auto generate3(const Eigen::Vector3d &C,
               const Eigen::Matrix<T, 3, Eigen::Dynamic> &xs)
    -> Eigen::Vector<T, Eigen::Dynamic> {
  std::vector<T> tmps;

  for (decltype(xs.cols()) i = 0; i < xs.cols(); ++i) {
    tmps.push_back(1.0 / (C - xs.col(i)).norm());
    // ys(i) = 1.0 / (C - xs.col(i)).norm();
  }

  Eigen::Vector<T, Eigen::Dynamic> ys =
      Eigen::Map<Eigen::Vector<T, Eigen::Dynamic>>(tmps.data(), tmps.size());

  return ys;
}
// The scalar function for which the gradient is needed
template <typename T>
T foo(const Eigen::Vector<T, Eigen::Dynamic> &x) {

  SPDLOG_INFO("x {}", x);
  return x.dot(x);
}

int main() {
  spdlog::set_pattern("[%Y-%m-%d %T][%^%l%$][TID: %t, PID: %P][%!@%s:%4#] %v");

  using T = dual;

  Eigen::Vector3d C(1.0, -1.0, M_PI);

  Eigen::Matrix<T, 3, Eigen::Dynamic> xs{{0.0, 0.0}, {0.0, 0.0}, {0.0, M_PI}};

  Eigen::Vector<T, 3> D(1.0, 0.0, -M_PI);

  //auto [f, dfdx] = derivatives(foo<T>, wrt(D), at(D));
  //SPDLOG_INFO("f {}", f);
  //SPDLOG_INFO("dfdx {}", dfdx);

  //auto [g, dfdxs] = derivatives(generate3<T>, wrt(xs), at(C, xs));
  //SPDLOG_INFO("dfdxs {}", dfdxs);

  Eigen::MatrixXd Z = jacobian(generate3<T>, wrt(xs), at(C, xs));
  SPDLOG_INFO("Z of dimension {}x{}", Z.rows(), Z.cols());
  SPDLOG_INFO("Z\n{}", Z);

  return EXIT_SUCCESS;
}
