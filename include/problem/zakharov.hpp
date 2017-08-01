#ifndef DE_OPTIMIZATION_PROBLEM_ZAKHAROV_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_ZAKHAROV_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class ZakharovFunction
 * \brief The Zakharov function
 *
 * \f$ f_{\mbox{Zakharov}}(\mathbf{x}) =
 * \displaystyle
 * \sum_{i=1}^{n} x_i^2 +
 * \left( \sum_{i=1}^D 0.5 i x_i \right)^2 +
 * \left( \sum_{i=1}^D 0.5 i x_i \right)^4 \f$
 */

class ZakharovFunction : public CECFunction<double> {
 public:
  explicit ZakharovFunction(const std::size_t D)
      : CECFunction(D, "Zakharov Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum_1 = 0.0, sum_2 = 0.0;
    for (std::size_t i = 0; i < D_; ++i) {
      sum_1 += pow(chromosome[i], 2);
      sum_2 += 0.5 * (i + 1) * chromosome[i];
    }
    return sum_1 + pow(sum_2, 2) + pow(sum_2, 4);
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_ZAKHAROV_FUNCTION_HPP
