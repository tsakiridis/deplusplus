#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_3_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_3_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/cigar.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/lunacek_bi_rastrigin.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction3
 * \brief Hybrid function 3
 */

class HybridFunction3 : public CECHybrid<double> {
 public:
  HybridFunction3(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 3",
                          {0.3, 0.3, 0.4},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {new CigarFunction(genes_[0]),
                  new RosenbrockFunction(genes_[1]),
                  new LunacekBiRastriginFunction(genes_[2], false)};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_3_HPP
