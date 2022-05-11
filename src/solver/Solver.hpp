#pragma once

#include <Eigen/Core>

namespace classyq {
/** Solver interface.
 *
 * @note We are using the curiously recurring template pattern (CRTP)
 */
template<typename T>
class Solver
{
public:
  Eigen::VectorXd compute(const Eigen::VectorXd& RHS, int irrep = 0) const
  {
    return static_cast<const T*>(this)->compute(RHS, irrep);
  }
};

/**
 * @tparam Derived
 * @tparam LHSOp
 * @tparam RHSOp
 */
template<typename Derived, typename LHSOp, typename RHSOp>
class Direct final
{
public:
  using value_type = typename Derived::PlainObject::Scalar;

private:
  Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> A_;
  Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> B_;

public:
  // at this point the matrices/matrix-vector product functions are already created and only need be applied
  Direct()
    : A_{ LHSOp() }
    , B_{ RHSOp() }
  {}
  Eigen::Vector<value_type, Eigen::Dynamic> compute(const Eigen::MatrixBase<Derived>& RHS, int irrep = 0) const
  {
    auto b = B_(RHS, irrep);
    return A_(b, irrep);
  }
};
} // namespace classyq
