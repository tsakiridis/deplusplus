/*!
 * This example demonstrates the basic use of DE++.
 *
 * First a problem and an algorithm are selected (here Griewank's function and
 * SHADE / L-SHADE respectively). The algorithms are ran, timed and the results
 * are presented in the standard output.
 */

#include "problem/griewank.hpp"
#include "algorithm/shade.hpp"
#include "timer.hpp"
#include <iostream>

size_t SEED = 100;

int main() {
  constexpr std::size_t D = 50;  // Dimension of the problem
  std::unique_ptr<DE::Problem::Base<double>> f;
  {
    std::cout << "50-dimensional Griewank function with SHADE" << std::endl;
    Timer t;
    f = std::make_unique<DE::Problem::GriewankFunction>(D);
    DE::Algorithm::SHADE<double> shade(std::move(f), false);  // SHADE
    shade.evolve_population();
    auto solution = shade.get_best();
    std::chrono::duration<double> elapsed_seconds = t.elapsed();
    std::cout << "Best fitness: " << solution.best_fitness << std::endl;
    std::cout << "Seconds elapsed: " << elapsed_seconds.count() << std::endl;
  }
  {
    std::cout << "50-dimensional Griewank function with L-SHADE" << std::endl;
    Timer t;
    f = std::make_unique<DE::Problem::GriewankFunction>(D);
    DE::Algorithm::SHADE<double> shade(std::move(f), true);  // L-SHADE
    shade.evolve_population();
    auto solution = shade.get_best();
    std::chrono::duration<double> elapsed_seconds = t.elapsed();
    std::cout << "Best fitness: " << solution.best_fitness << std::endl;
    std::cout << "Seconds elapsed: " << elapsed_seconds.count() << std::endl;
  }
}
