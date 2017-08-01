#ifndef DE_OPTIMIZATION_PROBLEM_SCHWEFEL_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_SCHWEFEL_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class SchwefelFunction
 * \brief Modified Schwefel's function
 *
 * \f$ \displaystyle
 * \begin{align}
 * f_{\mbox{Schwefel}}(\mathbf{x}) &=
 * 418.9829 D - \sum_{i=1}^{D} g(z_i) \\
 * z_i &= x_i + 4.209687462275036\mathrm{e}{+002} \\
 * g(z_i) &=
 * \begin{cases}
 * z_i \sin{\left ( \left | z_i \right |^{1/2} \right )}, &
 *  \text{if } |z_i| \leq 500 \\
 * (500 - \text{mod}(z_i, 500)
 * \sin{\left( \sqrt{|500 - \text{mod}(z_i, 500)|} \right)} -
 * \displaystyle \frac{(z_i - 500)^2}{10000 D}, & \text{if } z_i > 500 \\
 * (\text{mod}(z_i, 500 - 500)
 * \sin{\left( \sqrt{|\text{mod}(z_i, 500) - 500|} \right)} -
 * \displaystyle \frac{(z_i + 500)^2}{10000 D}, & \text{if } z_i < -500
 * \end{cases} \\
 * \end{align}
 * \f$
 */

class SchwefelFunction : public CECFunction<double> {
 public:
  explicit SchwefelFunction(const std::size_t D)
      : CECFunction(D, "Schwefel's Function") {
    scale_ = 1000.0 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0;
    for (std::size_t i = 0; i < D_; ++i) {
      double z = 4.209687462275036e+002 + chromosome[i];
      if (z > 500) {
        sum -= (500.0 - fmod(z, 500)) * sin(pow(500.0 - fmod(z, 500), 0.5)) -
               pow((z - 500.0) / 100, 2) / D_;
      } else if (z < -500) {
        sum -= (-500.0 + fmod(fabs(z), 500)) *
                   sin(pow(500.0 - fmod(fabs(z), 500), 0.5)) -
               pow((z + 500.0) / 100, 2) / D_;
      } else {
        sum -= z * sin(pow(fabs(z), 0.5));
      }
    }
    sum += 4.189828872724338e+002 * D_;
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_SCHWEFEL_FUNCTION_HPP
