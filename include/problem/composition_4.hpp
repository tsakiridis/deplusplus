#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_4_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_4_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/ackley.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/griewank.hpp"
#include "problem/rastrigin.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction4
 * \brief Composition function 4 of CEC-2017
 */

class CompositionFunction4 : public CECComposition<double> {
 public:
  CompositionFunction4(const std::size_t D,
                       const char* shift_file = nullptr,
                       const char* rotation_file = nullptr)
      : CECComposition<double>(D, "Composition Function 4") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0, 300.0};
    functions_ = {
        BasicFunction(new AckleyFunction(D), 10, 10, bias[0]),
        BasicFunction(new HighConditionedElliptic(D), 20, 1e-6, bias[1]),
        BasicFunction(new GriewankFunction(D), 30, 10, bias[2]),
        BasicFunction(new RastriginFunction(D), 40, 1, bias[3])};
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
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_4_HPP
