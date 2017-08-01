#ifndef DE_OPTIMIZATION_PROBLEM_HAPPYCAT_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_HAPPYCAT_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HappyCatFunction
 * \brief HappyCat function
 *
 * \f$ \displaystyle
 * f_{\mbox{HappyCat}}(\mathbf{x}) =
 * \left | \sum_{i=1}^{D} x_i^2 - D \right |^{1/4} + (0.5 \sum_{i=1}^{D} x_i^2 +
 * \sum_{i=1}^{D} x_i) / D + 0.5
 * \f$
 */

class HappyCatFunction : public CECFunction<double> {
 public:
  explicit HappyCatFunction(const std::size_t D)
      : CECFunction(D, "HappyCat Function") {
    scale_ = 5.0 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    /* original global optimum: [-1,-1,...,-1] */
    double squared_sum = 0.0, gene_sum = 0.0;
    constexpr double alpha = 1.0 / 8.0;
    for (std::size_t i = 0; i < D_; ++i) {
      squared_sum += (chromosome[i] - 1.0) * (chromosome[i] - 1.0);
      gene_sum += chromosome[i] - 1.0;
    }
    return pow(fabs(squared_sum - D_), 2 * alpha) +
           (0.5 * squared_sum + gene_sum) / D_ + 0.5;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HAPPYCAT_FUNCTION_HPP
