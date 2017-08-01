#ifndef DE_OPTIMIZATION_PROBLEM_GRIEWANK_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_GRIEWANK_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class GriewankFunction
 * \brief Griewank's function
 *
 * \f$ \displaystyle
 * f_{\mbox{Griewank}}(\mathbf{x}) =
 * 1 + \sum_{i=1}^{D} \frac{x_i^2}{400} -
 * \prod_{i=1}^{D} \cos{ \left( \frac{x_i}{\sqrt{i}} \right )}
 * \f$
 */

class GriewankFunction : public CECFunction<double> {
 public:
  explicit GriewankFunction(const std::size_t D)
      : CECFunction(D, "Griewank's Function") {
    scale_ = 600.0 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0, product = 1.0;
    for (std::size_t i = 0; i < D_; ++i) {
      sum += chromosome[i] * chromosome[i];
      product *= cos(chromosome[i] / sqrt(i + 1));
    }
    return 1 + sum / 4000 - product;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_GRIEWANK_FUNCTION_HPP
