#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_1_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_1_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/rastrigin.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction1
 * \brief Composition function 1 of CEC-2017
 */

class CompositionFunction1 : public CECComposition<double> {
 public:
  CompositionFunction1(const std::size_t D,
                       const char* shift_file = nullptr,
                       const char* rotation_file = nullptr)
      : CECComposition<double>(D, "Composition Function 1") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0};
    functions_ = {
        BasicFunction(new RosenbrockFunction(D), 10, 1, bias[0]),
        BasicFunction(new HighConditionedElliptic(D), 20, 1e-6, bias[1]),
        BasicFunction(new RastriginFunction(D), 30, 1, bias[2])};
    if (shift_file)
      for (std::size_t i = 0; i < functions_.size(); ++i)
        functions_[i].func->parse_shift_file(shift_file, i);
    if (rotation_file)
      for (std::size_t i = 0; i < functions_.size(); ++i)
        functions_[i].func->parse_rotation_file(rotation_file, i * D);
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_1_HPP
