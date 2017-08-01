/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Declarations of the DEGL algorithm.
 *
 * DEGL is a template class. Under the src/ folder you can find
 * the implementation for this class methods. Allowed templates
 * are float and double. Use double for extra precision.
 */

#ifndef DE_DEGL_HPP
#define DE_DEGL_HPP

#include <cstdint>
#include <vector>
#include <array>
#include "algorithm/base_algorithm.hpp"

namespace DE {
namespace Algorithm {

/*!
 *  \class DEGL
 *  \brief The DEGL algorithm
 *
 *  The DEGL algorithm \cite Das2009 (Das, S., Abraham, A., Chakraborty, U. K.,
 *  & Konar, A. (2009). Differential Evolution Using a Neighborhood-Based
 *  Mutation Operator. IEEE Transactions on Evolutionary Computation, 13(3))
 *  is a variant of the DE algorithm which can be used for the optimization of
 *  real-parameter problems. It utilizes the concept of small (local) neighbor-
 *  hoods for each population member, defined in a ring topology.
 *  To make the DEGL class work, a class must be defined inheriting the
 *  abstract base class Base, implementing the pure virtual functions of the
 *  class
 *
 *  \see Base
 *  \see SimpleFitnessFunction
 */

template <class T>
class DEGL : public Base<T> {
 public:
  /*!
   * \brief Create a new optimizer from an existing initial chromosome
   *
   * \param problem            : Pointer to a Base Problem
   * \param initial_chromosome : An initial (ideally not random) solution
   * \param minimize           : If true, minimize the fitness function
   */

  DEGL(std::shared_ptr<Problem::Base<T>> problem,
       const std::vector<T>& initial_chromosome,
       const bool minimize = true);

  /*!
   * \brief Create a new optimizer with no a priori knowledge
   *
   * \param problem  : Pointer to a Base Problem
   * \param minimize : If true, minimize the fitness function
   */

  DEGL(std::shared_ptr<Problem::Base<T>> problem, const bool minimize = true);

  /*!
   * \brief Apply DEGL
   *
   * \param max_generations : Maximum number of evolution generations
   *
   * \return The best solution
   */

  void evolve_population(const std::size_t max_generations = 5000);

 private:
  static constexpr float F = 0.8;  /*!< Scale factor */
  static constexpr float Cr = 0.9; /*!< Crossover factor */
  const std::size_t k_;            /*!< Neighborhood size */
  std::vector<float> w_;           /*!< Weight for every solution [N_] */
  std::vector<float> w_mutated_;   /*!< Mutated weight parameters [N_] */

  /*!
   * \brief Initialize the weights
   */

  void initialize_weights();

  /*!
   * Ring topology; find which indexes are within the neighborhood of a given
   * index.
   *
   * \param index  : The current position
   *
   * \return A vector with all the valid indexes
   */

  std::vector<std::size_t> get_neighborhood_indexes(
      const std::size_t index) const;

  /*!
   * Ring topology; get random indexes from within a neighborhood
   *
   * \param places  : The vector from get_valid_indexes
   *
   * \return An array with two random indexes.
   */

  std::array<std::size_t, 2> get_random_local_indexes(
      const std::vector<std::size_t>& places) const;

  /*!
   * \brief Get two random different indexes from [0, N_ - 1].
   *
   * Furthermore, they should be different from \param avoid
   *
   * \param avoid : the index to avoid
   */

  std::array<std::size_t, 2> get_random_global_indexes(
      const std::size_t avoid) const;

  /*!
   * Find the best from within a neighborhood.
   *
   * \param places : Vector with valid places (get_valid_indexes)
   *
   * \return The position of the local best
   */

  std::size_t find_local_best(const std::vector<std::size_t>& places) const;

  /*!
   * Mutation in DEGL [Das, et. al].
   *
   * Mutation as described in DEGL. Get 2 random local indexes and 2 global
   * indexes. Find the local best (global best is passed from outside).
   * Mutate w and chromosome and ensure constraints.
   *
   * \param index       : chromosome to be mutated
   * \param global_best : index of the global best chromosome
   */

  std::vector<T> mutate(const std::size_t index, const std::size_t global_best);
};  // class DEGL
}  // namespace Algorithm
}  // namespace DE
#endif  // DE_DEGL_HPP
