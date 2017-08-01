#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_5_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_5_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/cigar.hpp"
#include "problem/hgbat.hpp"
#include "problem/rastrigin.hpp"
#include "problem/rosenbrock.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction5
 * \brief Hybrid function 5
 */

class HybridFunction5 : public CECHybrid<double> {
 public:
  HybridFunction5(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 5",
                          {0.2, 0.2, 0.3, 0.3},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {new CigarFunction(genes_[0]), new HGBatFunction(genes_[1]),
                  new RastriginFunction(genes_[2]),
                  new RosenbrockFunction(genes_[3])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_5_HPP
