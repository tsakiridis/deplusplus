#ifndef DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_7_HPP
#define DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_7_HPP

#include "problem/cec_composition_function.hpp"
#include "problem/hgbat.hpp"
#include "problem/rastrigin.hpp"
#include "problem/modified_schwefel.hpp"
#include "problem/cigar.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/schaffer.hpp"

namespace DE {
namespace Problem {

/*!
 * \class CompositionFunction7
 * \brief Composition function 7 of CEC-2017
 */

class CompositionFunction7 : public CECComposition<double> {
 public:
  CompositionFunction7(const std::size_t D,
                       const char* shift_file = nullptr,
                       const char* rotation_file = nullptr)
      : CECComposition<double>(D, "Composition Function 7") {
    auto bias = std::vector<double>{0.0, 100.0, 200.0, 300.0, 400.0, 500.0};
    functions_ = {
        BasicFunction(new HGBatFunction(D), 10, 10, bias[0]),
        BasicFunction(new RastriginFunction(D), 20, 10, bias[1]),
        BasicFunction(new SchwefelFunction(D), 30, 2.5, bias[2]),
        BasicFunction(new CigarFunction(D), 40, 1e-26, bias[3]),
        BasicFunction(new HighConditionedElliptic(D), 50, 1e-6, bias[4]),
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
#endif  // DE_OPTIMIZATION_PROBLEM_COMPOSITION_FUNCTION_7_HPP
