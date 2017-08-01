#ifndef DE_OPTIMIZATION_PROBLEM_RASTRIGIN_NON_CONTIINUOUS_ROTATED_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_RASTRIGIN_NON_CONTIINUOUS_ROTATED_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class RastriginNonContinuousRotatedFunction
 * \brief Non-continuous Rotated Rastrigin's function
 *
 * \f$ \displaystyle
 * f_{\mbox{Non-continuous Rotated Rastrigin}}(\mathbf{x}) =
 * \sum_{i=1}^{D}\left( x_i^2 - 10 \cos(2 \pi x_i) + 10 \right)
 * \f$
 */

class RastriginNonContinuousRotatedFunction : public CECFunction<double> {
 public:
  explicit RastriginNonContinuousRotatedFunction(const std::size_t D)
      : CECFunction(D, "Non-continuous Rotated Rastrigin's Function") {
    scale_ = 5.12 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0;
    for (std::size_t i = 0; i < D_; ++i)
      sum += (chromosome[i] * chromosome[i] -
              10.0 * cos(2.0 * M_PI * chromosome[i]) + 10.0);
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_RASTRIGIN_NON_CONTIINUOUS_ROTATED_FUNCTION_HPP
