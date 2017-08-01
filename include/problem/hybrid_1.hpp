#ifndef DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_1_HPP
#define DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_1_HPP

#include "problem/cec_hybrid_function.hpp"
#include "problem/zakharov.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/rastrigin.hpp"

namespace DE {
namespace Problem {

/*!
 * \class HybridFunction1
 * \brief Hybrid function 1
 */

class HybridFunction1 : public CECHybrid<double> {
 public:
  HybridFunction1(const std::size_t D,
                  const char* shuffle_file = nullptr,
                  const std::size_t shuffle_offset = 0)
      : CECHybrid<double>(D,
                          "Hybrid Function 1",
                          {0.2, 0.4, 0.4},
                          shuffle_file,
                          shuffle_offset) {
    functions_ = {new ZakharovFunction(genes_[0]),
                  new RosenbrockFunction(genes_[1]),
                  new RastriginFunction(genes_[2])};
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_OPTIMIZATION_PROBLEM_HYBRID_FUNCTION_1_HPP
