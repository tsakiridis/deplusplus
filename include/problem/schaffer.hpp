#ifndef DE_OPTIMIZATION_PROBLEM_SCHAFFER_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_SCHAFFER_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class SchafferFunction
 * \brief Expanded Schaffer's function
 *
 * \f$ \displaystyle
 * \begin{align}
 * f_{\mbox{Schaffer}}(\mathbf{x}) &=
 * 0.5 + \displaystyle
 *  \frac{\sin\left(\sqrt{x^2+y^2}\right)^2 - 0.5}{1 + 0.001(x^2 + y^2)} \\
 * f_{\mbox{ExpandedSchaffer}}(\mathbf{x}) &=
 * f_{\mbox{Schaffer}}(x_1, x_2) +
 * f_{\mbox{Schaffer}}(x_2, x_3) +
 * \ldots
 * f_{\mbox{Schaffer}}(x_{D-}1, x_D) +
 * f_{\mbox{Schaffer}}(x_D, x_1)
 * \end{align}
 * \f$
 */

class SchafferFunction : public CECFunction<double> {
 public:
  explicit SchafferFunction(const std::size_t D)
      : CECFunction(D, "Schaffer's Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = 0.0, temp_1 = 0.0, temp_2 = 0.0;
    for (std::size_t i = 0; i < D_ - 1; ++i) {
      temp_1 = sin(sqrt(chromosome[i] * chromosome[i] +
                        chromosome[i + 1] * chromosome[i + 1]));
      temp_1 = temp_1 * temp_1;
      temp_2 = 1.0 +
               0.001 * (chromosome[i] * chromosome[i] +
                        chromosome[i + 1] * chromosome[i + 1]);
      sum += 0.5 + (temp_1 - 0.5) / (temp_2 * temp_2);
    }
    temp_1 = sin(sqrt(chromosome[D_ - 1] * chromosome[D_ - 1] +
                      chromosome[0] * chromosome[0]));
    temp_1 = temp_1 * temp_1;
    temp_2 = 1.0 +
             0.001 * (chromosome[D_ - 1] * chromosome[D_ - 1] +
                      chromosome[0] * chromosome[0]);
    sum += 0.5 + (temp_1 - 0.5) / (temp_2 * temp_2);
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_SCHAFFER_FUNCTION_HPP
