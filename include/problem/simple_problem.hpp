/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Definitions for simple function problems
 */

#ifndef DE_SIMPLE_PROBLEM_HPP
#define DE_SIMPLE_PROBLEM_HPP

#include <map>
#include <limits>
#include "problem/base_problem.hpp"
#include "rand.hpp"

namespace DE {
namespace Problem {
/*! \struct Constrain
 *  \brief A simple struct containing a [lower, upper] constrain
 */

template <typename T>
struct Constrain {
 public:
  T lower; /*!< Lower bound */
  T upper; /*!< Upper bound */
  Constrain(const T& l, const T& u) : lower(l), upper(u){};
};

/*! \class SimpleFitnessFunction
 *  \brief A wrapper of Base for simple fitness functions
 *
 *  If you don't need to do explicit work inside your randomize and
 *  constrain functions inherit this class, and populate \p constrains_
 *  in your constructor. Then only implement the fitness function.
 */

template <class T>
class SimpleFitnessFunction : public Base<T> {
 public:
  /*!
   * \brief Constructor for a simple function
   *
   * \param D    : The dimension of the problem
   * \param name : The function's name
   */

  SimpleFitnessFunction(const std::size_t D, const char* name)
      : Base<T>(D), name_(name){};

  void randomize(std::vector<T>& chromosome) const {
    assert(chromosome.size() == Base<T>::D_);
    for (std::size_t i = 0; i < Base<T>::D_; ++i) {
      const auto& c = SimpleFitnessFunction<T>::constrains_.find(i);
      if (c != SimpleFitnessFunction<T>::constrains_.end()) {  // exists
        chromosome[i] = rand_uniform_real(c->second.lower, c->second.upper);
      } else {
        chromosome[i] = rand_uniform_real(std::numeric_limits<T>::min(),
                                          std::numeric_limits<T>::max());
      }
    }
  }

  void constrain(std::vector<T>& chromosome) const {
    assert(chromosome.size() == Base<T>::D_);
    for (std::size_t i = 0; i < Base<T>::D_; ++i) {
      const auto& c = SimpleFitnessFunction<T>::constrains_.find(i);
      if (c != SimpleFitnessFunction<T>::constrains_.end()) {  // exists
        chromosome[i] =
            rand_within(chromosome[i], c->second.lower, c->second.upper);
        // If you want to clip the values, \see clip
      }
    }
  }

  /*!
   * \brief Get the functions name
   *
   * \return The \p name_
   */

  const char* get_name() const { return name_; }

 protected:
  const char* name_; /*!< Name of the fitness function */
  /*! Map indices (genes) to Constrain - if only one limit exists, replace
   * the other one with std::numeric_limit<T>::min() or max() */
  std::map<std::size_t, Constrain<T>> constrains_;

 private:
  /*!
   * \brief Clip a \p value between [\p lower, \p upper]
   *
   * \param value : Input value
   * \param lower : Left limit of the interval
   * \param upper : Right limit of the interval
   *
   * \return The new clipped value
   */

  T clip(const T& value, const T& lower, const T& upper) const {
    return std::max(lower, std::min(value, upper));
  }

  /*!
   * \brief If a \p value doesn't belong in the interval
   *        [\p lower, \p upper], randomize it within this interval
   *
   * \param value : Input value
   * \param lower : Left limit of the interval
   * \param upper : Right limit of the interval
   *
   * \return The new value within [lower, upper]
   */

  T rand_within(const T& value, const T& lower, const T& upper) const {
    return (value < lower || value > upper) ? rand_uniform_real(lower, upper)
                                            : value;
  }
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_SIMPLE_PROBLEM_HPP
