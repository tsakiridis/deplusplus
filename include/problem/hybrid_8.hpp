#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_8_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_8_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/ackley.hpp"
#include "problem/rastrigin.hpp"
#include "problem/hgbat.hpp"
#include "problem/discus.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction8
 * \brief Hybrid function 8
 */

class HybridFunction8 : public CECHybrid<double> {
 public:
  HybridFunction8(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 8",
                          {0.2, 0.2, 0.2, 0.2, 0.2},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {new HighConditionedElliptic(genes_[0]),
                  new AckleyFunction(genes_[1]),
                  new RastriginFunction(genes_[2]),
                  new HGBatFunction(genes_[3]), new DiscusFunction(genes_[4])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_8_HPP
