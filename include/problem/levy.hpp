#ifndef DE_OPTIMIZATION_PROBLEM_LEVY_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_LEVY_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class LevyFunction
 * \brief The Levy function
 *
 * \f$ f_{\mbox{Levy}}(\mathbf{x}) =
 * \displaystyle
 * \sin^2(\pi w_1) +
 * \left( \sum_{i=1}^{D-1} \left\{
 *    (w_i-1)^2 (1 + 10\sin^2(\pi w_i + 1)) \right\} \right) +
 * (w_D-1)^2 (1 + \sin^2(2\pi w_D)) \f$
 */

class LevyFunction : public CECFunction<double> {
 public:
  explicit LevyFunction(const std::size_t D)
      : CECFunction(D, "Levy Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    std::vector<double> w(D_, 0.0);
    for (size_t i = 0; i < D_; ++i)
      w[i] = 1.0 + (chromosome[i] - 1.0) / 4.0;

    double term1 = pow((sin(M_PI * w[0])), 2),
           term3 = pow((w[D_ - 1] - 1), 2) *
                   (1 + pow((sin(2 * M_PI * w[D_ - 1])), 2)),
           sum = 0.0;

    for (size_t i = 0; i < D_ - 1; ++i)
      sum += pow((w[i] - 1), 2) * (1 + 10 * pow((sin(M_PI * w[i] + 1)), 2));

    return term1 + sum + term3;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_LEVY_FUNCTION_HPP
