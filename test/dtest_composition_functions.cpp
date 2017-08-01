#include "gtest/gtest.h"
#include <memory>
#include "problem/cec_composition_function.hpp"
#include "problem/composition_1.hpp"
#include "problem/composition_2.hpp"
#include "problem/composition_3.hpp"
#include "problem/composition_4.hpp"
#include "problem/composition_5.hpp"
#include "problem/composition_6.hpp"
#include "problem/composition_7.hpp"
#include "problem/composition_8.hpp"
#include "problem/composition_9.hpp"
#include "problem/composition_10.hpp"
#include "test_utils.hpp"
#include "cec17_test_func.hpp"

void initialize_function(
    std::unique_ptr<DE::Problem::CECComposition<double>>& function,
    const std::size_t num,
    const std::size_t D,
    const bool use_rotation) {
  auto shift_str = shift_file(num), rotation_str = rotation_file(num, D),
       shuffle_str = shuffle_file(num, D);
  auto shift_f = shift_str.c_str(),
       rotation_f = (use_rotation ? rotation_str.c_str() : nullptr),
       shuffle_f = shuffle_str.c_str();
  switch (num) {
    case 21:
      function = std::make_unique<DE::Problem::CompositionFunction1>(
          D, shift_f, rotation_f);
      break;
    case 22:
      function = std::make_unique<DE::Problem::CompositionFunction2>(
          D, shift_f, rotation_f);
      break;
    case 23:
      function = std::make_unique<DE::Problem::CompositionFunction3>(
          D, shift_f, rotation_f);
      break;
    case 24:
      function = std::make_unique<DE::Problem::CompositionFunction4>(
          D, shift_f, rotation_f);
      break;
    case 25:
      function = std::make_unique<DE::Problem::CompositionFunction5>(
          D, shift_f, rotation_f);
      break;
    case 26:
      function = std::make_unique<DE::Problem::CompositionFunction6>(
          D, shift_f, rotation_f);
      break;
    case 27:
      function = std::make_unique<DE::Problem::CompositionFunction7>(
          D, shift_f, rotation_f);
      break;
    case 28:
      function = std::make_unique<DE::Problem::CompositionFunction8>(
          D, shift_f, rotation_f);
      break;
    case 29:
      function = std::make_unique<DE::Problem::CompositionFunction9>(
          D, shift_f, rotation_f, shuffle_f);
      break;
    case 30:
      function = std::make_unique<DE::Problem::CompositionFunction10>(
          D, shift_f, rotation_f, shuffle_f);
      break;
  }
}

namespace {

class CompositionFunctions : public ::testing::Test {
 protected:
  CompositionFunctions() { init_values(); }
};

TEST_F(CompositionFunctions, same_values_shifted_and_rotated) {
  std::unique_ptr<DE::Problem::CECComposition<double>> function;
  double fitness;
  for (std::size_t i = 21; i <= 30; ++i) {
    for (auto& x : x_tests) {
      initialize_function(function, i, x.size(), true);
      cec17_test_func(x.data(), &fitness, x.size(), 1, i);
      EXPECT_DOUBLE_EQ(function->fitness(x) + 100.0 * i, fitness)
          << function->get_name();
    }
  }
}

}  // namespace
