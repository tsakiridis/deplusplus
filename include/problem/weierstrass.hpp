#ifndef DE_OPTIMIZATION_PROBLEM_WEIERSTRASS_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_WEIERSTRASS_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class WeierstrassFunction
 * \brief Weierstrass's function
 *
 * \f$ \displaystyle
 * f_{\mbox{Weierstrass}}(\mathbf{x}) =
 * \sum_{i=1}^{D} \left ( \sum_{k=0}^{k_{max}} \left [ a^k
 * \cos{\left ( 2 \pi b^k (x_i + 0.5) \right )} \right ] \right ) -
 * D \sum_{k=0}^{k_{max}} \left [a^k \cos{\left (2 \pi b^k \cdot 0.5)}]
 * \f$
 */

class WeierstrassFunction : public CECFunction<double> {
 public:
  explicit WeierstrassFunction(const std::size_t D)
      : CECFunction(D, "Weierstrass's Function") {
    scale_ = 0.5 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum_2 = 0.0, total = 0.0;
    constexpr double a = 0.5, b = 3, k_max = 20;
    for (std::size_t i = 0; i < D_; ++i) {
      double sum_1 = 0.0;
      sum_2 = 0.0;
      for (std::size_t j = 0; j <= k_max; ++j) {
        sum_1 +=
            pow(a, j) * cos(2.0 * M_PI * pow(b, j) * (chromosome[i] + 0.5));
        sum_2 += pow(a, j) * cos(2.0 * M_PI * pow(b, j) * 0.5);
      }
      total += sum_1;
    }
    return total - D_ * sum_2;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_WEIERSTRASS_FUNCTION_HPP
