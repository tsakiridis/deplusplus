#ifndef DE_OPTIMIZATION_PROBLEM_EXPANDED_GRIEWANK_PLUS_ROSENBROCK_HPP
#define DE_OPTIMIZATION_PROBLEM_EXPANDED_GRIEWANK_PLUS_ROSENBROCK_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class ExpandedGriewankPlusRosenbrock
 * \brief Expanded Griewank's plus Rosenbrock's function
 *
 * \f$ \displaystyle
 * f_{\mbox{EGPR}}(\mathbf{x}) =
 * f_{\mbox{Griewank}}(f_{\mbox{Rosenbrock}}(x_1, x_2)) +
 * f_{\mbox{Griewank}}(f_{\mbox{Rosenbrock}}(x_2, x_3)) +
 * \ldots
 * f_{\mbox{Griewank}}(f_{\mbox{Rosenbrock}}(x_{D-1}, x_D)) +
 * f_{\mbox{Griewank}}(f_{\mbox{Rosenbrock}}(x_D, x_1))
 * \f$
 */

class ExpandedGriewankPlusRosenbrock : public CECFunction<double> {
 public:
  explicit ExpandedGriewankPlusRosenbrock(const std::size_t D)
      : CECFunction(D, "Expanded Griewank's plus Rosenbrokc's Function") {
    scale_ = 5.0 / 100.0;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double rosenbrock_first_sum = 0.0, rosenbrock_output = 0.0, sum = 0.0;
    for (std::size_t i = 0; i < D_ - 1; ++i) {
      rosenbrock_first_sum =
          (chromosome[i] + 1) * (chromosome[i] + 1) - (chromosome[i + 1] + 1);
      rosenbrock_output = 100.0 * rosenbrock_first_sum * rosenbrock_first_sum +
                          chromosome[i] * chromosome[i];
      sum += (rosenbrock_output * rosenbrock_output) / 4000.0 -
             cos(rosenbrock_output) + 1.0;
    }
    rosenbrock_first_sum = (chromosome[D_ - 1] + 1) * (chromosome[D_ - 1] + 1) -
                           (chromosome[0] + 1);
    rosenbrock_output = 100.0 * rosenbrock_first_sum * rosenbrock_first_sum +
                        chromosome[D_ - 1] * chromosome[D_ - 1];
    sum += (rosenbrock_output * rosenbrock_output) / 4000.0 -
           cos(rosenbrock_output) + 1.0;
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_EXPANDED_GRIEWANK_PLUS_ROSENBROCK_HPP
