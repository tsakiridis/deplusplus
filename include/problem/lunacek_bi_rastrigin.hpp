#ifndef DE_OPTIMIZATION_PROBLEM_LUNACEK_BI_RASTRIGIN_FUNCTION_HPP
#define DE_OPTIMIZATION_PROBLEM_LUNACEK_BI_RASTRIGIN_FUNCTION_HPP

#include <cmath>
#include <assert.h>
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*!
 * \class LunacekBiRastriginFunction
 * \brief The Lunacek bi-Rastrigin function
 *
 * \f$ f_{\mbox{Lunacek's bi-Rastrigin}}(\mathbf{x}) =
 * \displaystyle
 * \text{min}\Bigg[
 *   \sum^D_{i=1}\left(\hat{x}_i-\mu_0\right)^2 \ , \ 
 *   d \cdot D + s \cdot \sum^D_{i=1}\left(\hat{x}_i-\mu_1\right)^2\Bigg]
 * + 10 (D - \sum^D_{i=1} cos(2 \pi \hat{z}_i))
 * \\
 * \mu_0 = 2.5, \:
 * \mu_1 = -\sqrt{\frac{\displaystyle \mu^2_0-d}{\displaystyle s}}, \:
 * s = 1 - \frac{\displaystyle 1}{\displaystyle 2\sqrt{D+20}-8.2}, \:
 * d = 1
 * \\
 * y = \frac{\displaystyle 10 (\mathbf{x} - \mathbf{o})}{100}, \:
 * \hat{x_i} = 2 \text{sign}(x_i) y_i + \mu_0,\: \text{for } i=1,2,\ldots,D
 * \\
 * \hat{z} = \Lambda^{100}(\hat{\mathbf{x}} - \mu_0) 
 * \f$
 */

class LunacekBiRastriginFunction : public CECFunction<double> {
 public:
  LunacekBiRastriginFunction(const std::size_t D, const bool b = true)
      : CECFunction(D, "Lunacek bi-Rastrigin Function") {
    handle_shift_and_rotation_internally_ = b;
  }

  double evaluate(const std::vector<double>& chromosome) const {
    assert(chromosome.size() == D_);
    double mu_0 = 2.5, d = 1.0;
    double s = 1.0 - 1.0 / (2.0 * pow(D_ + 20.0, 0.5) - 8.2);
    double mu_1 = -pow((mu_0 * mu_0 - d) / s, 0.5);

    std::vector<double> shifted =
        handle_shift_and_rotation_internally_ ? shift(chromosome) : chromosome;

    for (std::size_t i = 0; i < D_; ++i)
      shifted[i] =
          2 * (chromosome[i] < 0 ? -1.0 : 1.0) * 10.0 / 100.0 * shifted[i];

    std::vector<double> zeta =
        handle_shift_and_rotation_internally_ ? rotate(shifted) : shifted;

    for (auto& gene : shifted)
      gene += mu_0;

    double sum_1 = 0.0, sum_2 = 0.0, sum_3 = 0.0;
    for (std::size_t i = 0; i < D_; ++i) {
      sum_1 += (shifted[i] - mu_0) * (shifted[i] - mu_0);
      sum_2 += (shifted[i] - mu_1) * (shifted[i] - mu_1);
    }

    for (std::size_t i = 0; i < D_; ++i)
      sum_3 += cos(2.0 * M_PI * zeta[i]);

    return std::min(sum_1, d * D_ + s * sum_2) + 10 * (D_ - sum_3);
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_LUNACEK_BI_RASTRIGIN_FUNCTION_HPP
