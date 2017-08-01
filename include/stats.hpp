/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Simple implementation for statistics of 1D vectors
 */

#ifndef DE_STATS_HPP
#define DE_STATS_HPP

#include <algorithm>
#include <numeric>
#include <vector>

/*!
 * \class Stats
 * \brief A class to extract statistics from 1D vectors
 */

template <typename T>
class Stats {
  /*! The data whose statistics will be extracted is a 1D vector */
  using data_type = std::vector<T>;

 public:
  /*!
   * \brief Explicit constructor using a 1D vector dataset
   *
   * \param d : the dataset
   */

  explicit Stats(const data_type& d) : data(d){};

  /*!
   * \brief Calculate the average value
   *
   * \return \f$\frac{1}{n}\sum_{i=1}^n x_{i,dim}\f$
   */

  double mean() const {
    return std::accumulate(data.cbegin(), data.cend(), 0.0) / data.size();
  }

  /*!
   * \brief Find the median value
   *
   * \return the median
   */

  T median() const {
    auto copy = data;
    std::sort(copy.begin(), copy.end());
    if (copy.size() % 2 == 0)
      return (copy[copy.size() / 2 - 1] + copy[copy.size() / 2]) / 2;
    else
      return copy[copy.size() / 2];
  }

  /*!
   * \brief Find the maximum value
   *
   * \return the maximum
   */

  T max() const { return *std::max_element(data.cbegin(), data.cend()); }

  /*!
   * \brief Find the minimum value
   *
   * \return the minimum
   */

  T min() const { return *std::min_element(data.cbegin(), data.cend()); }

  /*!
   * \brief Calculate the variance
   *
   * Uses Bessel's correction
   *
   * \return \f$\displaystyle \frac{\sum_{i=1}^n (x_{i,dim} - \mu)^2}{n-1} \f$
   */

  double variance() const {
    auto avg = mean();
    return std::accumulate(data.begin(), data.end(), 0.0,
                           [&avg](double& sum, const T& value) {
                             return sum + (value - avg) * (value - avg);
                           });
  }

  /*! \struct Statistics
   *  \brief An object containing all statistic measures
   */

  struct Statistics {
    double mean;     /*!< Average value */
    double median;   /*!< Median value */
    double variance; /*!< Variance */
    T min;           /*!< Minimum value */
    T max;           /*!< Maximum value */
  };

  /**
   * \brief Calculate all statistics for the vector
   *
   * \return A Statistics struct containing the results
   */

  Statistics calculate_all() const {
    Statistics s;
    s.mean = mean();
    s.median = median();
    s.variance = variance();
    s.min = min();
    s.max = max();
    return s;
  }

 private:
  const data_type data; /*!< data from which the stats will be calculated */
};

#endif  // DE_STATS_HPP
