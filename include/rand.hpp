/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Implementation of a random engine
 */

#ifndef RANDOM_ENGINE
#define RANDOM_ENGINE

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <assert.h>

extern size_t SEED; /*!< The seed must be defined externally */

/*!
 * \brief Generate a random float between (min, max) from a uniform distribution
 *
 * \param min : lower value
 * \param max : upper value
 *
 * \return A random double within (min, max)
 */

inline double rand_uniform_real(const double min, const double max) {
  assert(min < max);
  thread_local boost::random::mt19937 generator(SEED);
  boost::random::uniform_real_distribution<> dist(min, max);
  return dist(generator);
}

/**
 * \brief Generate a random int within [min, max] from a uniform distribution
 *
 * \param min : lower value
 * \param max : upper value
 *
 * \return A random integer within [min, max]
 */

inline int rand_uniform_int(const float min, const float max) {
  assert(min < max);
  thread_local boost::random::mt19937 generator(SEED);
  boost::random::uniform_int_distribution<> dist(min, max);
  return dist(generator);
}

/*!
 * \brief Generate a random float from a normal distribution
 *
 * \param mean  : mean value of the normal distribution
 * \param sigma : standard deviation of the normal distribution
 *
 * \return A random double
 */

inline double rand_normal(const double mean, const double sigma) {
  thread_local boost::random::mt19937 generator(SEED);
  boost::random::normal_distribution<> dist(mean, sigma);
  return dist(generator);
}

/*!
 * \brief Generate a random float from a Cauchy distribution
 *
 * \param mean  : mean value of the Cauchy distribution
 * \param sigma : standard deviation of the Cauchy distribution
 *
 * \return A random double
 */

inline double rand_cauchy(const double mean, const double sigma) {
  thread_local boost::random::mt19937 generator(SEED);
  boost::random::normal_distribution<> dist(mean, sigma);
  return dist(generator);
}

#endif  // RANDOM_ENGINE
