#ifndef DE_OPTIMIZATION_PROBLEM_DISCUS_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_DISCUS_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class DiscusFunction
 * \brief Discus function
 *
 * \f$ \displaystyle
 * f_{\mbox{discus}}(\mathbf{x}) =
 * 10^6 x_i^2 + \sum_{i=2}^{D} x_i^2
 * \f$
 */

class DiscusFunction : public CECFunction<double> {
 public:
  explicit DiscusFunction(const std::size_t D)
      : CECFunction(D, "Discus Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 1000000 * chromosome[0] * chromosome[0];
    for (std::size_t i = 1; i < D_; ++i)
      sum += chromosome[i] * chromosome[i];
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_DISCUS_FUNCTION_HPP
