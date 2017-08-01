#include "gtest/gtest.h"
#include <memory>
#include "problem/cec_basic_problem.hpp"
#include "problem/cigar.hpp"
#include "problem/sum_of_different_power.hpp"
#include "problem/zakharov.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/rastrigin.hpp"
#include "problem/schaffer_f7.hpp"
#include "problem/lunacek_bi_rastrigin.hpp"
#include "problem/rastrigin_non_continuous_rotated.hpp"
#include "problem/levy.hpp"
#include "problem/modified_schwefel.hpp"
#include "test_utils.hpp"
#include "cec17_test_func.hpp"

void initialize_unimodal_function(
    std::unique_ptr<DE::Problem::CECFunction<double>>& function,
    const std::size_t num,
    const std::size_t D) {
  switch (num) {
    case 1:
      function.reset(new DE::Problem::CigarFunction(D));
      break;
    case 2:
      function.reset(new DE::Problem::SumOfDifferentPowerFunction(D));
      break;
    case 3:
      function.reset(new DE::Problem::ZakharovFunction(D));
      break;
    case 4:
      function.reset(new DE::Problem::RosenbrockFunction(D));
      break;
    case 5:
      function.reset(new DE::Problem::RastriginFunction(D));
      break;
    case 6:
      function.reset(new DE::Problem::SchafferF7Function(D));
      break;
    case 7:
      function.reset(new DE::Problem::LunacekBiRastriginFunction(D));
      break;
    case 8:
      function.reset(new DE::Problem::RastriginNonContinuousRotatedFunction(D));
      break;
    case 9:
      function.reset(new DE::Problem::LevyFunction(D));
      break;
    case 10:
      function.reset(new DE::Problem::SchwefelFunction(D));
      break;
  }
  function->parse_shift_file(shift_file(num).c_str());
  function->parse_rotation_file(rotation_file(num, D).c_str());
}

namespace {

class UnimodalFunctions : public ::testing::Test {
 protected:
  UnimodalFunctions() { init_values(); }
};

TEST_F(UnimodalFunctions, same_values) {
  std::unique_ptr<DE::Problem::CECFunction<double>> function;
  double fitness;
  for (std::size_t i = 1; i <= 10; ++i) {
    for (auto& x : x_tests) {
      initialize_unimodal_function(function, i, x.size());
      cec17_test_func(x.data(), &fitness, x.size(), 1, i);
      EXPECT_DOUBLE_EQ(function->fitness(x) + 100.0 * i, fitness)
          << function->get_name();
    }
  }
}

}  // namespace
