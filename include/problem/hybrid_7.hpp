#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_7_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_7_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/katsuura.hpp"
#include "problem/ackley.hpp"
#include "problem/expanded_griewank_plus_rosenbrock.hpp"
#include "problem/modified_schwefel.hpp"
#include "problem/rastrigin.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction7
 * \brief Hybrid function 7
 */

class HybridFunction7 : public CECHybrid<double> {
 public:
  HybridFunction7(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 7",
                          {0.1, 0.2, 0.2, 0.2, 0.3},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {
        new KatsuuraFunction(genes_[0]), new AckleyFunction(genes_[1]),
        new ExpandedGriewankPlusRosenbrock(genes_[2]),
        new SchwefelFunction(genes_[3]), new RastriginFunction(genes_[4])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_7_HPP
