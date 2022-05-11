#include "SolverImpl.hpp"

#include <Eigen/Core>

#include "cavity/TsLess.hpp"

namespace classyq {
[[nodiscard]] auto
form_S(const TsLess& cavity) -> Eigen::MatrixXd
{
  auto N = cavity.size();

  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N, N);

  auto kstart = 0;

  for (auto A = 0; A < cavity.n_spheres(); ++A) {
    auto sph_A = cavity.spheres(A);

    // total number of points in EQ partition of sphere A
    auto N_A = sph_A.size();
    // radius of sphere A
    auto R_A = sph_A.radius();
    // EQ weight of unit sphere derived from sphere A
    auto omega_A = sph_A.omega();
    // EQ points of unit sphere derived from sphere A
    const auto EQs = sph_A.points();

    // compute self-potential factors for unit sphere derived from sphere A
    Eigen::VectorXd fs = Eigen::VectorXd::Constant(N_A, 4.0 * M_PI);
    for (auto i = 0; i < N_A; ++i) {
      auto p_i = EQs.col(i);
      for (auto j = 0; j < N_A; ++j) {
        auto p_j = EQs.col(j);
        if (i != j) {
          auto dist = (p_i - p_j).norm();
          fs(i) -= (omega_A) / dist;
        }
      }
    }

    // number of points in the *final* quadrature belonging to sphere A
    auto m_A = cavity.points_per_sphere(A);
    S.diagonal()(Eigen::seqN(kstart, m_A)) = R_A * fs(Eigen::seqN(0, m_A));
    kstart += m_A;
  }

  const auto points = cavity.points();

  // fill the off-diagonal of S
  for (auto i = 0; i < N; ++i) {
    const auto s_i = points.col(i);
    for (auto j = i + 1; j < N; ++j) {
      const auto s_j = points.col(j);
      S(i, j) = S(j, i) = 1.0 / (s_i - s_j).norm();
    }
  }

  return S;
}

[[nodiscard]] auto
form_D(const TsLess& cavity) -> Eigen::MatrixXd
{
  auto N = cavity.size();

  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(N, N);

  auto kstart = 0;

  for (auto A = 0; A < cavity.n_spheres(); ++A) {
    auto sph_A = cavity.spheres(A);

    // total number of points in EQ partition of sphere A
    auto N_A = sph_A.size();
    // EQ weight of unit sphere derived from sphere A
    auto omega_A = sph_A.omega();
    // EQ points of unit sphere derived from sphere A
    const auto EQs = sph_A.points();

    // compute self-field factors for unit sphere derived from sphere A
    Eigen::VectorXd gs = Eigen::VectorXd::Constant(N_A, -2.0 * M_PI);
    for (auto i = 0; i < N_A; ++i) {
      auto p_i = EQs.col(i);
      for (auto j = 0; j < N_A; ++j) {
        auto p_j = EQs.col(j);
        if (i != j) {
          auto dist = (p_i - p_j).norm();
          gs(i) -= (omega_A * (p_i - p_j).dot(p_j)) / std::pow(dist, 3);
        }
      }
    }

    // number of points in the *final* quadrature belonging to sphere A
    auto m_A = cavity.points_per_sphere(A);
    D.diagonal()(Eigen::seqN(kstart, m_A)) = gs(Eigen::seqN(0, m_A));
    kstart += m_A;
  }

  const auto points = cavity.points();
  const auto normals = cavity.normals();

  // fill the off-diagonal of D
  for (auto i = 0; i < N; ++i) {
    auto p_i = points.col(i);
    for (auto j = 0; j < N; ++j) {
      auto p_j = points.col(j);
      auto n_j = normals.col(j);
      if (i != j) {
        auto dist = (p_i - p_j).norm();
        D(i, j) = (p_i - p_j).dot(n_j) / std::pow(dist, 3);
      }
    }
  }

  return D;
}

[[nodiscard]] auto
form_R_epsilon(const TsLess& cavity, double epsilon) -> Eigen::MatrixXd
{
  auto N = cavity.size();

  auto D = form_D(cavity);

  auto W = cavity.weights().asDiagonal();

  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(N, N);

  auto f_eps = (epsilon + 1.0) / (epsilon - 1.0);

  return (2 * M_PI * f_eps * Id - D * W);
}

[[nodiscard]] auto
form_R_infinity(const TsLess& cavity) -> Eigen::MatrixXd
{
  auto N = cavity.size();

  auto D = form_D(cavity);

  auto W = cavity.weights().asDiagonal();

  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(N, N);

  return (2 * M_PI * Id - D * W);
}
} // namespace classyq
