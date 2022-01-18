#pragma once

#include <Eigen/Core>

template <typename DerivedA, typename DerivedB>
bool allclose(const Eigen::DenseBase<DerivedA> &a,
              const Eigen::DenseBase<DerivedB> &b,
              const typename DerivedA::RealScalar &rtol = 1e-5,
              const typename DerivedA::RealScalar &atol = 1e-8) {
  return ((a.derived() - b.derived()).array().abs() <=
          (atol + rtol * b.derived().array().abs()))
      .all();
}
