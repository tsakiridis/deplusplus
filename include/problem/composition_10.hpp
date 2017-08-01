#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_10_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_10_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/hybrid_5.hpp"
#include "problem/hybrid_8.hpp"
#include "problem/hybrid_9.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction10
 * \brief Composition function 10 of CEC-2017
 */

class CompositionFunction10 : public CECComposition<double> {
 public:
  CompositionFunction10(const std::size_t D,
                        const char* shift_file = nullptr,
                        const char* rotation_file = nullptr,
                        const char* shuffle_file = nullptr)
      : CECComposition<double>(D, "Composition Function 10") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0};
    functions_ = {
        BasicFunction(new HybridFunction5(D, shuffle_file), 10, 1, bias[0]),
        BasicFunction(new HybridFunction8(D, shuffle_file, D), 30, 1, bias[1]),
        BasicFunction(new HybridFunction9(D, shuffle_file, 2 * D), 50, 1,
                      bias[2])};
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
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_10_HPP
