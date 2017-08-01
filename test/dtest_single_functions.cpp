#include "gtest/gtest.h"
#include <memory>
#include "problem/cec_basic_problem.hpp"
#include "problem/cigar.hpp"
#include "problem/sum_of_different_power.hpp"
#include "problem/zakharov.hpp"
#include "problem/rosenbrock.hpp"
#include "problem/rastrigin.hpp"
#include "problem/schaffer.hpp"
#include "problem/levy.hpp"
#include "problem/modified_schwefel.hpp"
#include "problem/high_conditioned_elliptic_function.hpp"
#include "problem/discus.hpp"
#include "problem/ackley.hpp"
#include "problem/weierstrass.hpp"
#include "problem/griewank.hpp"
#include "problem/katsuura.hpp"
#include "problem/happycat.hpp"
#include "problem/hgbat.hpp"
#include "problem/expanded_griewank_plus_rosenbrock.hpp"
#include "problem/schaffer_f7.hpp"
#include "test_utils.hpp"
#include "cec17_test_func.hpp"

void initialize_single_function(
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
      function.reset(new DE::Problem::SchafferFunction(D));
      break;
    // cases 7, 8 have built-in shifts and rotations
    case 9:
      function.reset(new DE::Problem::LevyFunction(D));
      break;
    case 10:
      function.reset(new DE::Problem::SchwefelFunction(D));
      break;
    case 11:
      function.reset(new DE::Problem::HighConditionedElliptic(D));
      break;
    case 12:
      function.reset(new DE::Problem::DiscusFunction(D));
      break;
    case 13:
      function.reset(new DE::Problem::AckleyFunction(D));
      break;
    case 14:
      function.reset(new DE::Problem::WeierstrassFunction(D));
      break;
    case 15:
      function.reset(new DE::Problem::GriewankFunction(D));
      break;
    case 16:
      function.reset(new DE::Problem::KatsuuraFunction(D));
      break;
    case 17:
      function.reset(new DE::Problem::HappyCatFunction(D));
      break;
    case 18:
      function.reset(new DE::Problem::HGBatFunction(D));
      break;
    case 19:
      function.reset(new DE::Problem::ExpandedGriewankPlusRosenbrock(D));
      break;
    case 20:
      function.reset(new DE::Problem::SchafferF7Function(D));
      break;
  }
}

void calculate_cec_single_function(double* x,
                                   double* fitness,
                                   std::size_t D,
                                   const std::size_t function) {
  switch (function) {
    case 1:
      bent_cigar_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 2:
      sum_diff_pow_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 3:
      zakharov_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 4:
      rosenbrock_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 5:
      rastrigin_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 6:
      escaffer6_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    // cases 7, 8 have built-in shifts and rotations
    case 9:
      levy_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 10:
      schwefel_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 11:
      ellips_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 12:
      discus_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 13:
      ackley_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 14:
      weierstrass_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 15:
      griewank_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 16:
      katsuura_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 17:
      happycat_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 18:
      hgbat_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 19:
      grie_rosen_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
    case 20:
      schaffer_F7_func(x, fitness, D, nullptr, nullptr, 0, 0);
      break;
  }
}

namespace {

class SimpleFunctions : public ::testing::Test {
 protected:
  SimpleFunctions() { init_values(); }

  static const std::vector<std::size_t> valid_indices;
  static const std::vector<std::size_t> zero_indices;
};

const std::vector<std::size_t> SimpleFunctions::zero_indices = {
    1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15};
const std::vector<std::size_t> SimpleFunctions::valid_indices = {
    1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15};

TEST_F(SimpleFunctions, zero_value) {
  std::unique_ptr<DE::Problem::CECFunction<double>> function;
  const std::vector<double> x{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (const auto& i : zero_indices) {
    initialize_single_function(function, i, x.size());
    EXPECT_NEAR(function->fitness(x), 0, 1e-14) << function->get_name();
  }
}

TEST_F(SimpleFunctions, Out_of_bounds) {
  std::unique_ptr<DE::Problem::CECFunction<double>> function;
  for (const auto& i : valid_indices) {
    initialize_single_function(function, i, 150);
    ASSERT_DEATH(function->evaluate({0, 0}), ".*chromosome.size\\(\\) == D_.*");
    ASSERT_DEATH(function->evaluate({0, 0, 0, 0, 0, 0}),
                 ".*chromosome.size\\(\\) == D_.*");
  }
}

TEST_F(SimpleFunctions, same_values) {
  std::unique_ptr<DE::Problem::CECFunction<double>> function;
  double fitness;
  for (const auto& i : valid_indices) {
    for (auto& x : x_tests) {
      initialize_single_function(function, i, x.size());
      calculate_cec_single_function(x.data(), &fitness, x.size(), i);
      EXPECT_DOUBLE_EQ(function->fitness(x), fitness) << function->get_name();
    }
  }
}

}  // namespace
