/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Definitions for CEC composition functions
 */

#ifndef DE_CEC_COMPOSITION_PROBLEM_HPP
#define DE_CEC_COMPOSITION_PROBLEM_HPP

#include <vector>
#include <numeric>
#include <assert.h>
#include <cstdio>
#include <limits>
#include <algorithm>
#include "problem/simple_problem.hpp"
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

template <class T>
class CECComposition : public CECFunction<T> {
 public:
  CECComposition(const std::size_t D, const char* name)
      : CECFunction<T>(D, name) {
    // Constrains
    for (std::size_t i = 0; i < Base<T>::D_; ++i)
      SimpleFitnessFunction<T>::constrains_.insert(
          std::make_pair(i, Constrain<double>(-100, 100)));
  };

  double evaluate(const std::vector<T>& chromosome) const {
    assert(chromosome.size() == Base<T>::D_);
    thread_local std::vector<double> fitness, weights;
    fitness.resize(functions_.size());
    weights.resize(functions_.size());
    for (std::size_t i = 0; i < functions_.size(); ++i) {
      fitness[i] =
          functions_[i].lambda * functions_[i].func->fitness(chromosome) +
          functions_[i].bias;
      const auto& shift = functions_[i].func->get_shift_data();
      assert(shift.size() == Base<T>::D_);
      weights[i] = 0;
      for (std::size_t j = 0; j < Base<T>::D_; ++j)
        weights[i] += pow(chromosome[j] - shift[j], 2.0);
      if (weights[i] == 0)
        weights[i] = std::numeric_limits<double>::max();
      else
        weights[i] = pow(1.0 / weights[i], 0.5) *
                     exp(-1.0 * weights[i] / 2.0 / Base<T>::D_ /
                         pow(functions_[i].sigma, 2.0));
    }
    double weights_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (weights_sum == 0.0) {
      for (auto& w : weights)
        w = 1.0;
      weights_sum = functions_.size();
    }
    double fit_sum = 0.0;
    for (std::size_t i = 0; i < functions_.size(); ++i)
      fit_sum += weights[i] / weights_sum * fitness[i];
    return fit_sum;
  }

 protected:
  struct BasicFunction {
    CECFunction<double>* func;
    double sigma;
    double lambda;
    double bias;

    BasicFunction(CECFunction<double>* f, double s, double l, double b)
        : func(f), sigma(s), lambda(l), bias(b){};
  };

  std::vector<BasicFunction> functions_;
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_CEC_COMPOSITION_PROBLEM_HPP
