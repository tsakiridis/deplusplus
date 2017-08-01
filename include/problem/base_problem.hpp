/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Contains the basic interface for all problems
 */

#ifndef DE_BASE_PROBLEM_HPP
#define DE_BASE_PROBLEM_HPP

#include <vector>

namespace DE {

/*!
 * \namespace Problem
 * \brief All optimization problems should be inside this namespace
 */

namespace Problem {

/*! \class Base
 *  \brief The problem to be optimized
 *
 *  To use any of the DE algorithms, an optimization problem must be
 * defined. Create a class inheriting this abstract base class and
 * implement the pure virtual functions.
 *
 *   For simple fitness functions a wrapper is provided which assists
 *  you in the process. \see SimpleFitnessFunction
 */

template <class T>
class Base {
 public:
  /*!
   * \brief Simple constructor
   *
   * Every problem has a dimensionality!
   *
   * \param D : Dimension of chromosome, number of genes
   */

  explicit Base(const std::size_t D) : D_(D){};

  /*!
   * \brief Randomize a given chromosome
   *
   * \param chromosome : the chromosome to be randomized
   */

  virtual void randomize(std::vector<T>& chromosome) const = 0;

  /*!
   * \brief Constrain a given chromosome
   *
   * Respect the constrains for each different gene
   *
   * \param chromosome : the chromosome to be constrained
   */

  virtual void constrain(std::vector<T>& chromosome) const = 0;

  /*!
   * \brief Calculate the fitness of the chromosome
   *
   * \param chromosome : the chromosome to be evaluated
   */

  virtual double fitness(const std::vector<T>& chromosome) const = 0;

  /*!
   * \brief Get the number of genes
   *
   * \return The dimensionality of the problem
   */

  std::size_t get_number_of_genes() const { return D_; }

 protected:
  const std::size_t D_; /*!< Number of genes */
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_BASE_PROBLEM_HPP
