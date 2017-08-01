#ifndef DE_OPTIMIZATION_PROBLEM_HIGH_CONDITIONED_ELLIPTIC_HPP
#define DE_OPTIMIZATION_PROBLEM_HIGH_CONDITIONED_ELLIPTIC_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HighConditionedElliptic
 * \brief High Conditioned Elliptic function
 *
 * \f$ \displaystyle
 * f_{\mbox{elliptic}}(\mathbf{x}) =
 * \sum_{i=1}^{D} (10^6) ^{\frac{i-1}{D-1}} x_i^2
 * \f$
 */

class HighConditionedElliptic : public CECFunction<double> {
 public:
  explicit HighConditionedElliptic(const std::size_t D)
      : CECFunction(D, "High Conditioned Elliptic Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0;
    for (std::size_t i = 0; i < D_; ++i)
      sum += pow(10, 6.0 * i / (D_ - 1)) * chromosome[i] * chromosome[i];
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HIGH_CONDITIONED_ELLIPTIC_HPP
