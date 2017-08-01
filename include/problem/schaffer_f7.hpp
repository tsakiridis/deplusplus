#ifndef DE_OPTIMIZATION_PROBLEM_SCHAFFER_F7_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_SCHAFFER_F7_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class SchafferF7Function
 * \brief Schaffer's F7 function
 *
 * \f$ \displaystyle
 * f_{\mbox{Schaffer-F7}}(\mathbf{x}) =
 * \Bigg[\frac{1}{D-1}\sum_{i=1}^{D-1}
 *    \sqrt{\displaystyle s_i} \cdot (sin(50.0s_i^{0.2})+1)\Bigg]^2, \:
 * s_i = \sqrt{x_i^2 + x_{i+1}^2}
 * \f$
 */

class SchafferF7Function : public CECFunction<double> {
 public:
  explicit SchafferF7Function(const std::size_t D)
      : CECFunction(D, "Schaffer's F7 Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0;
    for (std::size_t i = 0; i < D_ - 1; ++i) {
      double s = pow(
          chromosome[i] * chromosome[i] + chromosome[i + 1] * chromosome[i + 1],
          0.5);
      double tmp = sin(50.0 * pow(s, 0.2));
      sum += pow(s, 0.5) + pow(s, 0.5) * tmp * tmp;
    }
    return sum * sum / (D_ - 1) / (D_ - 1);
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_SCHAFFER_F7_FUNCTION_HPP
