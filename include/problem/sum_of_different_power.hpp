#ifndef DE_OPTIMIZATION_PROBLEM_SUM_OF_DIFFERENT_POWER_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_SUM_OF_DIFFERENT_POWER_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class SumOfDifferentPowerFunction
 * \brief Sum of different power function
 *
 * \f$ \displaystyle
 * f_{\mbox{sum_of_different_power}}(\mathbf{x}) =
 * \sum_{i=1}^{D} |x_i|^{i+1}
 * \f$
 */

class SumOfDifferentPowerFunction : public CECFunction<double> {
 public:
  explicit SumOfDifferentPowerFunction(const std::size_t D)
      : CECFunction(D, "Sum of different power Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0;
    for (std::size_t i = 0; i < D_; ++i)
      sum += pow((std::abs(chromosome[i])), (i + 1));
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_SUM_OF_DIFFERENT_POWER_FUNCTION_HPP
