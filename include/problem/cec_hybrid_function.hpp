/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Definitions for CEC hybrid function problems
 */

#ifndef DE_CEC_HYBRID_PROBLEM_HPP
#define DE_CEC_HYBRID_PROBLEM_HPP

#include <vector>
#include <numeric>
#include <assert.h>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include "problem/simple_problem.hpp"
#include "problem/cec_basic_problem.hpp"

namespace DE {
namespace Problem {

/*! \class CECHybrid
 *  \brief Class for defining CEC Hybrid functions
 */

template <class T>
class CECHybrid : public CECFunction<T> {
 public:
  CECHybrid(const std::size_t D,
            const char* name,
            const std::vector<double>& percentage,
            const char* shuffle_file = nullptr,
            const std::size_t shuffle_offset = 0)
      : CECFunction<T>(D, name), percentage_(percentage) {
    // Constrains
    for (std::size_t i = 0; i < Base<T>::D_; ++i)
      SimpleFitnessFunction<T>::constrains_.insert(
          std::make_pair(i, Constrain<double>(-100, 100)));
    // Populate genes_
    genes_.resize(percentage_.size(), 0);
    for (std::size_t i = 0; i < percentage_.size() - 1; ++i)
      genes_[i] = ceil(percentage_[i] * Base<T>::D_);
    genes_.back() =
        Base<T>::D_ - std::accumulate(genes_.begin(), genes_.end(), 0);
    parse_shuffle_data(shuffle_file, shuffle_offset);
  };

  double evaluate(const std::vector<T>& chromosome) const {
    assert(chromosome.size() == Base<T>::D_);
    thread_local std::vector<std::vector<T>> partials;
    partials.resize(percentage_.size());
    for (std::size_t i = 0; i < percentage_.size(); ++i)
      partials[i].resize(genes_[i]);
    double sum = 0.0;
    std::size_t start = 0;
    for (std::size_t i = 0; i < percentage_.size(); ++i) {
      start += (i == 0 ? 0 : genes_[i - 1]);
      const std::size_t end =
          (i == percentage_.size() ? Base<T>::D_ : start + genes_[i]);
      for (std::size_t j = start; j < end; ++j) {
        partials[i][j - start] = chromosome[shuffle_[j] - 1];
      }
      sum += functions_[i]->fitness(partials[i]);
    }
    return sum;
  }

  /*!
   * \brief Parse the shuffle file (shuffle_*.txt) generated by the .m file
   *
   * \param shuffle_file   : Path to the .txt file
   * \param shuffle_offset : If specified, skip this many values
   *
   * \throw std::runtime_error if the file cannot be opened
   */

  void parse_shuffle_data(const char* shuffle_file,
                          std::size_t shuffle_offset = 0) {
    shuffle_.resize(Base<T>::D_);
    if (shuffle_file) {
      FILE* fpt = fopen(shuffle_file, "r");
      if (fpt == NULL)
        throw std::runtime_error("Unable to open shuffle file");
      if (shuffle_offset > 0)
        while (shuffle_offset-- != 0)
          fscanf(fpt, "%lu", &shuffle_[0]);
      for (std::size_t i = 0; i < Base<T>::D_; ++i)
        fscanf(fpt, "%lu", &shuffle_[i]);
      fclose(fpt);
    } else {
      std::iota(shuffle_.begin(), shuffle_.end(), Base<T>::D_);
    }
  }

 protected:
  /*! The basic functions from which the hybrid function is comprised */
  std::vector<CECFunction<double>*> functions_;
  /*! The number of genes per function */
  std::vector<unsigned short int> genes_;

 private:
  /*! Percentage for each basic function */
  const std::vector<double> percentage_;
  /*! Random shuffling for all the genes (size D) */
  std::vector<std::size_t> shuffle_;
};
}  // namespace Problem
}  // namespace DE
#endif  // DE_CEC_HYBRID_PROBLEM_HPP
