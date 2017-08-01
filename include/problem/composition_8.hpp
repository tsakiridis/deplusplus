#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_8_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_8_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/ackley.hpp"
#include "problem/griewank.hpp"
#include "problem/discus.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/happycat.hpp"
#include "problem/schaffer.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction8
 * \brief Composition function 8 of CEC-2017
 */

class CompositionFunction8 : public CECComposition<double> {
 public:
  CompositionFunction8(const std::size_t D,
                       const char* shift_file = nullptr,
                       const char* rotation_file = nullptr)
      : CECComposition<double>(D, "Composition Function 8") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0, 300.0, 400.0, 500.0};
    functions_ = {BasicFunction(new AckleyFunction(D), 10, 10, bias[0]),
                  BasicFunction(new GriewankFunction(D), 20, 10, bias[1]),
                  BasicFunction(new DiscusFunction(D), 30, 1e-6, bias[2]),
                  BasicFunction(new RosenbrockFunction(D), 40, 1, bias[3]),
                  BasicFunction(new HappyCatFunction(D), 50, 1, bias[4]),
                  BasicFunction(new SchafferFunction(D), 60, 5e-4, bias[5])};
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
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_8_HPP
