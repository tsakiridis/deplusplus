#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_6_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_6_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/schaffer.hpp"
#include "problem/hgbat.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/modified_schwefel.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction6
 * \brief Hybrid function 6
 */

class HybridFunction6 : public CECHybrid<double> {
 public:
  HybridFunction6(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 6",
                          {0.2, 0.2, 0.3, 0.3},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {new SchafferFunction(genes_[0]), new HGBatFunction(genes_[1]),
                  new RosenbrockFunction(genes_[2]),
                  new SchwefelFunction(genes_[3])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_6_HPP
