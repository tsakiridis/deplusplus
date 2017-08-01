#ifndef DE_OPTIMIZATION_PROBLEM_KATSUURA_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_KATSUURA_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class KatsuuraFunction
 * \brief Katsuura's function
 *
 * \f$ \displaystyle
 * f_{\mbox{Katsuura}}(\mathbf{x}) =
 * \frac{10}{D^2} \prod_{i=1}^{D} \left(1 + i \sum_{j=1}^{32} \frac{|2^j x_i -
 * \text{round}(2^j x_i)|}{2^j} \right )^{\displaystyle \frac{10}{D^{1.2}}} -
 * \frac{10}{D^2}
 * \f$
 */

class KatsuuraFunction : public CECFunction<double> {
 public:
  explicit KatsuuraFunction(const std::size_t D)
      : CECFunction(D, "Katsuura Function") {
    scale_ = 5.0 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double product = 1.0, exponent = 10.0 / pow(D_, 1.2);
    for (std::size_t i = 0; i < D_; ++i) {
      double sum = 0.0;
      for (std::size_t j = 1; j <= 32; ++j) {
        double power_of_two = pow(2.0, j),
               inner_product = power_of_two * chromosome[i];
        sum += fabs(inner_product - round(inner_product)) / power_of_two;
      }
      product *= pow(1.0 + (i + 1) * sum, exponent);
    }
    double temp = 10.0 / D_ / D_;
    return temp * product - temp;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_KATSUURA_FUNCTION_HPP
