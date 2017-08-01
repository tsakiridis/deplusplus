/*!
 * \file
 * \author Nikos Tsakiridis <tsakirin@auth.gr>
 * \version 1.0
 *
 * \brief Simple implementation for a timer function
 */

#include <chrono>

/*!
 * \struct Timer
 * \brief A simple wrapper around chrono to record the time
 */

struct Timer {
  /*! The clock to be used */
  using clock = std::chrono::high_resolution_clock;

  /*! Count time elapsed */
  template <class DurationType = std::chrono::nanoseconds>
  auto elapsed() const {
    return std::chrono::duration_cast<DurationType>(clock::now() - start_);
  }

  /*! Reset the timer */
  void reset() { start_ = clock::now(); }

 private:
  /*! Clock starts ticking at construction time */
  std::chrono::time_point<clock> start_ = clock::now();
};
