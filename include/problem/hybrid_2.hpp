#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_2_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_2_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/modified_schwefel.hpp"
#include "problem/cigar.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction2
 * \brief Hybrid function 2
 */

class HybridFunction2 : public CECHybrid<double> {
 public:
  HybridFunction2(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 2",
                          {0.3, 0.3, 0.4},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {new HighConditionedElliptic(genes_[0]),
                  new SchwefelFunction(genes_[1]),
                  new CigarFunction(genes_[2])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_2_HPP
