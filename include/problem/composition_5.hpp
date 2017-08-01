#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_5_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_5_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/rastrigin.hpp"
#include "problem/happycat.hpp"
#include "problem/ackley.hpp"
#include "problem/discus.hpp"
#include "problem/rosenbrock.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction5
 * \brief Composition function 5 of CEC-2017
 */

class CompositionFunction5 : public CECComposition<double> {
 public:
  CompositionFunction5(const std::size_t D,
                       const char* shift_file = nullptr,
                       const char* rotation_file = nullptr)
      : CECComposition<double>(D, "Composition Function 5") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0, 300.0, 400.0};
    functions_ = {BasicFunction(new RastriginFunction(D), 10, 10, bias[0]),
                  BasicFunction(new HappyCatFunction(D), 20, 1, bias[1]),
                  BasicFunction(new AckleyFunction(D), 30, 10, bias[2]),
                  BasicFunction(new DiscusFunction(D), 40, 1e-6, bias[3]),
                  BasicFunction(new RosenbrockFunction(D), 50, 1, bias[4])};
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
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_5_HPP
