#ifndef DE_OPTIMIZATION_PROBLEM_CIGAR_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_CIGAR_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CigarFunction
 * \brief The Cigar function
 *
 * \f$ f_{\mbox{Cigar}}(\mathbf{x}) =
 * x_1^2 + \displaystyle 10^6 \sum_{i=2}^{n} x_i^2 \f$
 */

class CigarFunction : public CECFunction<double> {
 public:
  explicit CigarFunction(const std::size_t D)
      : CECFunction(D, "Cigar Function") {}

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double sum = chromosome[0] * chromosome[0];
    for (std::size_t i = 1; i < D_; ++i)
      sum += 1000000 * chromosome[i] * chromosome[i];
    return sum;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_CIGAR_FUNCTION_HPP
