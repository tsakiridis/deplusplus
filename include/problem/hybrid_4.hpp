#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_4_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_4_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/ackley.hpp"
#include "problem/schaffer_f7.hpp"
#include "problem/rastrigin.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction4
 * \brief Hybrid function 4
 */

class HybridFunction4 : public CECHybrid<double> {
 public:
  HybridFunction4(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 4",
                          {0.2, 0.2, 0.2, 0.4},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {
        new HighConditionedElliptic(genes_[0]), new AckleyFunction(genes_[1]),
        new SchafferF7Function(genes_[2]), new RastriginFunction(genes_[3])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_4_HPP
