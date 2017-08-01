#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_10_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_10_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/hgbat.hpp"
#include "problem/katsuura.hpp"
#include "problem/ackley.hpp"
#include "problem/rastrigin.hpp"
#include "problem/modified_schwefel.hpp"
#include "problem/schaffer_f7.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction10
 * \brief Hybrid function 10
 */

class HybridFunction10 : public CECHybrid<double> {
 public:
  HybridFunction10(const std::size_t D,
                   const char* shuffle_file = nullptr,
                   const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 10",
                          {0.1, 0.1, 0.2, 0.2, 0.2, 0.2},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {
        new HGBatFunction(genes_[0]),    new KatsuuraFunction(genes_[1]),
        new AckleyFunction(genes_[2]),   new RastriginFunction(genes_[3]),
        new SchwefelFunction(genes_[4]), new SchafferF7Function(genes_[5])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_10_HPP
