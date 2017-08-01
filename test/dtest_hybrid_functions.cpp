#include "gtest/gtest.h"
#include <memory>
#include "problem/cec_hybrid_function.hpp"
#include "problem/hybrid_1.hpp"
#include "problem/hybrid_2.hpp"
#include "problem/hybrid_3.hpp"
#include "problem/hybrid_4.hpp"
#include "problem/hybrid_5.hpp"
#include "problem/hybrid_6.hpp"
#include "problem/hybrid_7.hpp"
#include "problem/hybrid_8.hpp"
#include "problem/hybrid_9.hpp"
#include "problem/hybrid_10.hpp"
#include "test_utils.hpp"
#include "cec17_test_func.hpp"

void initialize_function(
    std::unique_ptr<DE::Problem::CECHybrid<double>>& function,
    const std::size_t num,
    const std::size_t D) {
  switch (num) {
    case 11:
      function = std::make_unique<DE::Problem::HybridFunction1>(
          D, shuffle_file(num, D).c_str());
      break;
    case 12:
      function = std::make_unique<DE::Problem::HybridFunction2>(
          D, shuffle_file(num, D).c_str());
      break;
    case 13:
      function = std::make_unique<DE::Problem::HybridFunction3>(
          D, shuffle_file(num, D).c_str());
      break;
    case 14:
      function = std::make_unique<DE::Problem::HybridFunction4>(
          D, shuffle_file(num, D).c_str());
      break;
    case 15:
      function = std::make_unique<DE::Problem::HybridFunction5>(
          D, shuffle_file(num, D).c_str());
      break;
    case 16:
      function = std::make_unique<DE::Problem::HybridFunction6>(
          D, shuffle_file(num, D).c_str());
      break;
    case 17:
      function = std::make_unique<DE::Problem::HybridFunction7>(
          D, shuffle_file(num, D).c_str());
      break;
    case 18:
      function = std::make_unique<DE::Problem::HybridFunction8>(
          D, shuffle_file(num, D).c_str());
      break;
    case 19:
      function = std::make_unique<DE::Problem::HybridFunction9>(
          D, shuffle_file(num, D).c_str());
      break;
    case 20:
      function = std::make_unique<DE::Problem::HybridFunction10>(
          D, shuffle_file(num, D).c_str());
      break;
  }
}

namespace {

class HybridFunctions : public ::testing::Test {
 protected:
  HybridFunctions() { init_values(); }
};

TEST_F(HybridFunctions, same_values_shifted) {
  std::unique_ptr<DE::Problem::CECHybrid<double>> function;
  double fitness;
  for (std::size_t i = 11; i <= 20; ++i) {
    for (auto& x : x_tests) {
      initialize_function(function, i, x.size());
      function->parse_shift_file(shift_file(i).c_str());
      function->parse_rotation_file(rotation_file(i, x.size()).c_str());
      cec17_test_func(x.data(), &fitness, x.size(), 1, i);
      EXPECT_DOUBLE_EQ(function->fitness(x) + 100.0 * i, fitness)
          << function->get_name();
    }
  }
}

}  // namespace
