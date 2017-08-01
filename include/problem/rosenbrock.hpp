#ifndef DE_OPTIMIZATION_PROBLEM_ROSENBROCK_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_ROSENBROCK_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class RosenbrockFunction
 * \brief Rosenbrock's function
 *
 * \f$ \displaystyle
 * f_{\mbox{Rosenbrock}}(\mathbf{x}) =
 * \sum_{i=1}^{D-1}(100(x_i^2-x_{i+1})^2 + (x_i - 1)^2)
 * \f$
 */

class RosenbrockFunction : public CECFunction<double> {
 public:
  explicit RosenbrockFunction(const std::size_t D)
      : CECFunction(D, "Rosenbrock's Function") {
    scale_ = 2.048 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0;
    for (std::size_t i = 0; i < D_ - 1; ++i)
      sum += 100 * pow((chromosome[i] + 1.0) * (chromosome[i] + 1.0) -
                           (chromosome[i + 1] + 1.0),
                       2) +
             pow(chromosome[i], 2);
    // Slightly modified (every gene + 1.0) to shift the minimum at 0
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_ROSENBROCK_FUNCTION_HPP
