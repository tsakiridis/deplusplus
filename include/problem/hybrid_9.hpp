#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_9_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_9_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/cigar.hpp"
#include "problem/rastrigin.hpp"
#include "problem/expanded_griewank_plus_rosenbrock.hpp"
#include "problem/weierstrass.hpp"
#include "problem/schaffer.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction9
 * \brief Hybrid function 9
 */

class HybridFunction9 : public CECHybrid<double> {
 public:
  HybridFunction9(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 9",
                          {0.2, 0.2, 0.2, 0.2, 0.2},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {
        new CigarFunction(genes_[0]), new RastriginFunction(genes_[1]),
        new ExpandedGriewankPlusRosenbrock(genes_[2]),
        new WeierstrassFunction(genes_[3]), new SchafferFunction(genes_[4])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_9_HPP
