/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Implementation of the SHADE 1.1 and L-SHADE algorithms.
 *
 * SHADE \cite Tanabe2013 is an adaptive DE algorithm which incorporates
 * success-history based parameter adaptation and one of the state-of-the-art
 * DE algorithms. L-SHADE \cite Tanabe2014 further extends
 * SHADE with Linear Population Size Reduction (LPSR), which
 * continually decreases the population size according to a linear
 * function.
 */

#ifndef L_SHADE_HPP
#define L_SHADE_HPP

/*! Special value set in the crossover memory values */
#define TERMINAL_VALUE -100

#include <vector>
#include <mutex>
#include "algorithm/base_algorithm.hpp"

namespace DE {
namespace Algorithm {

/*!
 * \brief The core SHADE implementation
 *
 * To make the SHADE class work, a class must be defined inheriting the
 * abstract base class Base, implementing the pure virtual functions of the
 * class
 *
 * \see Base
 * \see SimpleFitnessFunction
 */

template <class T>
class SHADE : public Base<T> {
 public:
  /*!
   * \brief Create a new optimizer from an existing initial chromosome
   *
   * \param problem                   : Pointer to a Base Problem
   * \param initial_chromosome        : An initial (ideally not random) solution
   * \param use_linear_size_reduction : If true, use L-SHADE
   * \param minimize                  : If true, minimize the fitness function
   */

  SHADE(std::shared_ptr<Problem::Base<T>> problem,
        const std::vector<T>& initial_chromosome,
        const bool use_linear_size_reduction = false,
        const bool minimize = true);

  /*!
   * \brief Create a new optimizer from an existing initial chromosome and its
   *        fitness
   *
   * \param problem                   : Pointer to a Base Problem
   * \param initial_chromosome        : An initial (ideally not random) solution
   * \param initial_fitness           : Fitness of initial_chromosome
   * \param use_linear_size_reduction : If true, use L-SHADE
   * \param minimize                  : If true, minimize the fitness function
   */

  SHADE(std::shared_ptr<Problem::Base<T>> problem,
        const std::vector<T>& initial_chromosome,
        const double initial_fitness,
        const bool use_linear_size_reduction = false,
        const bool minimize = true);

  /*!
   * \brief Create a new optimizer with no a priori knowledge
   *
   * \param problem                   : Pointer to a Base Problem
   * \param use_linear_size_reduction : If true, use L-SHADE
   * \param minimize                  : If true, minimize the fitness function
   */

  SHADE(std::shared_ptr<Problem::Base<T>> problem,
        const bool use_linear_size_reduction = false,
        const bool minimize = true);

  /*!
   * \brief Apply the SHADE algorithm
   *
   * \param max_generations : Maximum number of generations
   *
   * \return The best population
   */

  void evolve_population(const std::size_t max_generations = 5000);

 private:
  std::size_t N_;                        /*!< Number of chromosomes */
  std::size_t p_;                        /*!< p in current-to-pbest */
  const std::size_t H_;                  /*!< Size of memory */
  std::size_t A_size_;                   /*!< Maximum size of archive */
  std::vector<std::vector<T>> A_;        /*!< Chromosome archive */
  std::vector<std::size_t> top_p_;       /*!< Indices of the top p solutions */
  std::vector<float> Cr_;                /*!< Crossover memory values (H_) */
  std::vector<float> F_;                 /*!< Scale factor memory values (H_) */
  const bool use_linear_size_reduction_; /*!< If true, use L-SHADE */
  std::mutex pop_mutex_;                 /*!< Mutex for parallel processing */

  /*!
   * \brief Generate the crossover factor
   *
   * Refer to Equation (1) in \cite Tanabe2014
   *
   * \param rand_index : Index between [0, H_) pointing to Cr_
   *
   * \return The Cr to be used
   */

  float get_crossover_factor(const std::size_t rand_index) const;

  /*!
   * \brief Generate the scale factor
   *
   * Refer to Equation (2) in \cite Tanabe2014
   *
   * \param rand_index : Index between [0, H_) pointing to F_
   *
   * \return The F to be used
   */

  float get_scale_factor(const std::size_t rand_index) const;

  /*!
   * \brief Update the memory values for Cr and F
   *
   * Refer to Algorithm 1 in \cite Tanabe2014
   *
   * \param S_Cr      : Set of successful Crossover values
   * \param S_F       : Set of successful Scale factor values
   * \param delta_fit : delta f for each successful value
   * \param k         : Index to the memory place to be updated
   */

  void memory_update(const std::vector<float>& S_Cr,
                     const std::vector<float>& S_F,
                     const std::vector<double>& delta_fit,
                     const std::size_t k);

  /*!
   * \brief Calculate the weighted Lehmer mean
   *
   * \param weights : Weights for each value
   * \param values  : The values
   *
   * \return \f$\frac{\displaystyle\sum_{k} w_k \times (v_k)^2}
   *              {\displaystyle\sum_{k} w_k \times v_k}\f$
   */

  float weighted_lehmer_mean(const std::vector<double>& weights,
                             const std::vector<float>& values) const;

  /*!
   * \brief Add a \p chromosome to the archive, but respect A_size_
   *
   * \param chromosome : The chromosome to be appended to A_
   */

  void add_to_archive(const std::vector<T>& chromosome);

  /*!
   * \brief Find the top p solutions in x_ and update top_p_
   */

  void update_top_p_solutions();

  /*!
   * \brief Mutate a chromosome using current-to-pbest/1
   *
   * \param base_index : The current index
   * \param F          : The scale factor
   *
   * \return The mutant
   */

  std::vector<T> mutate(const std::size_t base_index, const float F) const;

  /*!
   * \brief Apply linear size reduction
   *
   * Instead of using FE as in the paper, we use the generation index
   *
   * \param current_generation : Current generation index
   * \param max_generations    : The max number of generations
   */

  void linear_size_reduction(const std::size_t current_generation,
                             const std::size_t max_generations);

  void parallel_evolve(const std::size_t thread_id,
                       std::vector<float>& S_Cr,
                       std::vector<float>& S_F,
                       std::vector<double>& delta_fit);

};  // class SHADE
}  // namespace Algorithm
}  // namespace DE

#endif  // L_SHADE_HPP
