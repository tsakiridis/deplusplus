/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Includes all CEC-2017 functions
 */

#include <memory>
#include <string>
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

/*! Base path for all the .txt files, set externally */
extern std::string base_path;

/*!
 * \brief Get the path to the shift file
 *
 * \param problem : Function number
 */

std::string shift_file(const std::size_t problem) {
  return base_path + "shift_data_" + std::to_string(problem) + ".txt";
}

/*!
 * \brief Get the path to the rotation file
 *
 * \param problem : Function number
 * \param D       : Dimensions
 */

std::string rotation_file(const std::size_t problem, const std::size_t D) {
  return base_path + "M_" + std::to_string(problem) + "_D" + std::to_string(D) +
         ".txt";
}

/*!
 * \brief Get the path to the shuffle file
 *
 * \param problem : Function number
 * \param D       : Dimensions
 */

std::string shuffle_file(const std::size_t problem, const std::size_t D) {
  return base_path + "shuffle_data_" + std::to_string(problem) + "_D" +
         std::to_string(D) + ".txt";
}

/*!
 * \brief Create a new function for CEC-2017
 *
 * \param function : Pointer to the function to be created
 * \param num      : Which function to create [1-15]
 * \param D        : Dimensions of the problem
 */

void initialize_function(
    std::unique_ptr<DE::Problem::CECFunction<double>>& function,
    const std::size_t num,
    const std::size_t D) {
  auto shift_str = shift_file(num), rotation_str = rotation_file(num, D),
       shuffle_str = shuffle_file(num, D);
  auto shift_f = shift_str.c_str(), rotation_f = rotation_str.c_str(),
       shuffle_f = shuffle_str.c_str();
  switch (num) {
    case 1:
      function = std::make_unique<DE::Problem::CigarFunction>(D);
      break;
    case 2:
      function = std::make_unique<DE::Problem::SumOfDifferentPowerFunction>(D);
      break;
    case 3:
      function = std::make_unique<DE::Problem::ZakharovFunction>(D);
      break;
    case 4:
      function = std::make_unique<DE::Problem::RosenbrockFunction>(D);
      break;
    case 5:
      function = std::make_unique<DE::Problem::RastriginFunction>(D);
      break;
    case 6:
      function = std::make_unique<DE::Problem::SchafferFunction>(D);
      break;
    case 7:
      function = std::make_unique<DE::Problem::LunacekBiRastriginFunction>(D);
      break;
    case 8:
      function =
          std::make_unique<DE::Problem::RastriginNonContinuousRotatedFunction>(
              D);
      break;
    case 9:
      function = std::make_unique<DE::Problem::LevyFunction>(D);
      break;
    case 10:
      function = std::make_unique<DE::Problem::SchwefelFunction>(D);
      break;
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
  if (num >= 1 && num <= 10) {  // Uni- and multi-modal
    function->parse_shift_file(shift_file(num).c_str());
    function->parse_rotation_file(rotation_file(num, D).c_str());
  }
}
