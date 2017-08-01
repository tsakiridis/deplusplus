#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_2_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_2_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/rastrigin.hpp"
#include "problem/griewank.hpp"
#include "problem/modified_schwefel.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction2
 * \brief Composition function 2 of CEC-2017
 */

class CompositionFunction2 : public CECComposition<double> {
 public:
  CompositionFunction2(const std::size_t D,
                       const char* shift_file = nullptr,
                       const char* rotation_file = nullptr)
      : CECComposition<double>(D, "Composition Function 2") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0};
    functions_ = {BasicFunction(new RastriginFunction(D), 10, 1, bias[0]),
                  BasicFunction(new GriewankFunction(D), 20, 10, bias[1]),
                  BasicFunction(new SchwefelFunction(D), 30, 1, bias[2])};
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
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_2_HPP
