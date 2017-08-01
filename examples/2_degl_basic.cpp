/*!
 * This example demonstrates the basic use of DE++.
 *
 * First a problem and an algorithm are selected; we select here the first
 * Hyrid function from CEC-2017 without the use of shuffling. The algorithm
 * tested is DEGL. The function is optimized, and the results are displayed
 * in standard output.
 */

#include "problem/hybrid_1.hpp"
#include "algorithm/degl.hpp"
#include "timer.hpp"
#include <iostream>

size_t SEED = 100;

int main() {
  constexpr std::size_t D = 10;
  std::unique_ptr<DE::Problem::Base<double>> f;
  Timer t;
  f = std::make_unique<DE::Problem::HybridFunction1>(D);
  DE::Algorithm::DEGL<double> degl(std::move(f), true);
  degl.evolve_population();
  auto solution = degl.get_best();
  std::chrono::duration<double> elapsed_seconds = t.elapsed();
  std::cout << solution.best_fitness << std::endl;
  std::cout << "Seconds elapsed: " << elapsed_seconds.count() << std::endl;
}
