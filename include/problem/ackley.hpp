#ifndef DE_OPTIMIZATION_PROBLEM_ACKLEY_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_ACKLEY_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class AckleyFunction
 * \brief Ackley's function
 *
 * \f$ \displaystyle
 * f_{\mbox{Ackley}}(\mathbf{x}) =
 * -20 \exp\left(-0.2 \sqrt{\frac{1}{D} \sum_{i=1}^{D} x_i^2}\right)
 * -\exp\left(\frac{1}{D} \sum_{i=1}^{D} \cos(2 \pi x_i)\right) + 20 + \exp(1)
 * \f$
 */

class AckleyFunction : public CECFunction<double> {
 public:
  explicit AckleyFunction(const std::size_t D)
      : CECFunction(D, "Ackley's Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum_1 = 0.0, sum_2 = 0.0;
    for (std::size_t i = 0; i < D_; ++i) {
      sum_1 += chromosome[i] * chromosome[i];
      sum_2 += cos(2 * M_PI * chromosome[i]);
    }
    return -20.0 * std::exp(-0.2 * sqrt(sum_1 / D_)) - std::exp(sum_2 / D_) +
           20 + std::exp(1.0);
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_ACKLEY_FUNCTION_HPP
