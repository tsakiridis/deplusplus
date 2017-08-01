/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Implementation of an abstract DE algorithm.
 */

#ifndef DE_BASE_ALGORITHM_HPP
#define DE_BASE_ALGORITHM_HPP

#include <cstddef>
#include <algorithm>
#include <assert.h>
#include <memory>
#include "problem/base_problem.hpp"

namespace DE {

/*!
 * \namespace Algorithm
 * \brief All DE algorithms should be inside this namespace
 */

namespace Algorithm {

/*!
 * \class Base
 * \brief Implementation of an abstract DE algorithm
 *
 *  DE algorithms should inherit this class.
 */

template <typename T>
class Base {
 public:
  /*!
   * \brief Create a new optimizer with no a priori knowledge
   *
   * \param problem  : Pointer to a Base Problem
   * \param N        : Number of chromosomes
   * \param minimize : If true, minimize the fitness function
   */

  Base(std::shared_ptr<Problem::Base<T>> problem,
       const std::size_t N,
       const bool minimize)
      : p_problem_(problem),
        D_(p_problem_->get_number_of_genes()),
        N_(N),
        minimize_(minimize),
        allow_parallel_(false) {
    fit_.resize(N_);
    std::vector<T> initial_chromosome(D_, 0);
    x_.resize(N_, initial_chromosome);
    for (std::size_t i = 0; i < N_; ++i) {
      p_problem_->randomize(x_[i]);
      fit_[i] = p_problem_->fitness(x_[i]);
    }
  };

  /*!
   * \brief Create a new optimizer from an existing initial chromosome
   *
   * \param problem            : Pointer to a Base Problem
   * \param initial_chromosome : An initial (ideally not random) solution
   * \param N                  : Number of chromosomes
   * \param minimize           : If true, minimize the fitness function
   */

  Base(std::shared_ptr<Problem::Base<T>> problem,
       const std::vector<T>& initial_chromosome,
       const std::size_t N,
       const bool minimize)
      : p_problem_(problem),
        D_(p_problem_->get_number_of_genes()),
        N_(N),
        minimize_(minimize),
        allow_parallel_(false) {
    assert(initial_chromosome.size() == D_);
    x_.resize(N_, initial_chromosome);
    fit_.resize(N_);
    fit_[0] = p_problem_->fitness(x_[0]);
    for (std::size_t i = 1; i < N_; ++i) {
      p_problem_->randomize(x_[i]);
      fit_[i] = p_problem_->fitness(x_[i]);
    }
  };

  /*!
   * \brief Main function body for every derived class
   *
   * \param max_generations : Maximum number of generations
   */

  virtual void evolve_population(const std::size_t max_generations) = 0;

  /*! \struct Results
   *  \brief A struct containing the best individual
   */

  struct Results {
    std::vector<T> best_chromosome; /*!< Best result */
    double best_fitness;            /*!< Best fitness */
    Results(const std::vector<T>& c, const double f)
        : best_chromosome(c), best_fitness(f){};
  };

  /*!
   * \brief Find and return the best result
   *
   * \return A Results object with the best individual
   */

  Results get_best() const {
    const auto best = best_index();
    return Results(x_[best], fit_[best]);
  }

  /*!
   * \brief Enables parallel (multi-thread) processing
   *
   * Call this function after the constructor to enable parallel computations
   */

  void allow_parallel_computations() { allow_parallel_ = true; }

 protected:
  /*! Class containing the fitness function to be optimized */
  const std::shared_ptr<Problem::Base<T>> p_problem_;
  const std::size_t D_;           /*!< Genes (dimensions of each chromosome) */
  const std::size_t N_;           /*!< Number of initial chromosomes */
  const bool minimize_;           /*!< If true, minimize the fitness function */
  std::vector<std::vector<T>> x_; /*!< Chromosome population (N_, D_) */
  std::vector<double> fit_;       /*!< Fitness of each solution (N_) */
  bool allow_parallel_;           /*!< True to enable parallel computations */

  /*!
   * \brief Ascertain if rhs fitness isn't worse than lhs
   *
   * \param lhs : The old fitness
   * \param rhs : The new fitness
   *
   * \return True if rhs isn't worse than lhs
   */

  bool compare_fitnesses_with_equality(const double& lhs,
                                       const double& rhs) const {
    return minimize_ ? lhs >= rhs : lhs <= rhs;
  }

  /*!
   * \brief Ascertain if rhs fitness is better than lhs
   *
   * \param lhs : The old fitness
   * \param rhs : The new fitness
   *
   * \return True if rhs is better than lhs
   */

  bool compare_fitnesses(const double& lhs, const double& rhs) const {
    return minimize_ ? lhs > rhs : lhs < rhs;
  }

  /*!
   * \brief Find the best individual in the population
   *
   * \return Index of the best individual
   */

  std::size_t best_index() const {
    return std::distance(
        fit_.begin(), (minimize_ ? std::min_element(fit_.begin(), fit_.end())
                                 : std::max_element(fit_.begin(), fit_.end())));
  }
};

}  // namespace Algorithm
}  // namespace DE

#endif  // DE_BASE_ALGORITHM_HPP
