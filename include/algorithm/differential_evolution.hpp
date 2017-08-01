/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Useful global functions for various DE algorithms
 *
 * The core of the Differential Evolution (\cite Storn1997) algorithm is the
 * same for most of the variants.
 */

#ifndef DE_ENGINE_HPP
#define DE_ENGINE_HPP

#include <vector>
#include "rand.hpp"

namespace DE {
namespace Algorithm {

/*!
 * \brief Perform binary crossover
 *
 * \param target : The target (inital) chromosome
 * \param donor  : The donor (mutated) chromosome
 * \param Cr     : The crossover factor
 *
 * \return The trial vector
 */

template <typename T>
std::vector<T> binary_crossover(const std::vector<T>& target,
                                const std::vector<T>& donor,
                                const float& Cr) {
  auto trial = target;
  const std::size_t D = target.size(), j_rand = rand_uniform_int(0, D - 1);
  for (std::size_t j = 0; j < D; ++j)
    if (j == j_rand || rand_uniform_real(0, 1) <= Cr)
      trial[j] = donor[j];
  return trial;
}

/*!
 * \brief Perform exponential crossover
 *
 * \param target : The target (inital) chromosome
 * \param donor  : The donor (mutated) chromosome
 * \param Cr     : The crossover factor
 *
 * \return The trial vector
 */

template <typename T>
std::vector<T> exponential_crossover(const std::vector<T>& target,
                                     const std::vector<T>& donor,
                                     const float& Cr) {
  auto trial = target;
  const std::size_t D = target.size(), start = rand_uniform_int(0, D - 1);
  std::size_t L = 0;  // Length of the two point crossover
  do {
    ++L;
  } while (rand_uniform_real(0, 1) < Cr && L < D);
  for (std::size_t j = start; j < start + L; ++j)
    trial[start % D] = donor[start % D];
  return trial;
}
}  // namespace Algorithm
}  // namespace DE
#endif  // DE_ENGINE_HPP
