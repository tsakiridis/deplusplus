/*!
 * This example demonstrates the use of DE++ for CEC benchmark functions.
 *
 * If you haven't already, you need to generate the shift, rotation etc. files
 * using the Matlab / Octave script. Kindly check README.md.
 *
 * To define the base path where these files are located you need to edit the
 * global variable base_path, defined below.
 *
 * This example uses L-SHADE and runs all experiments; it further accumulates
 * the results and writes them out to a csv file.
 *
 * Please note that this process is very, *very*, slow as it runs
 * 15 experiments * 4 Dimensions * 51 times = 3060 experiments
 */

#include "problem/cec_all_functions.hpp"
#include "algorithm/shade.hpp"
#include "stats.hpp"
#include <fstream>
#include <iostream>
#include <string>

size_t SEED = 100;
std::string base_path = "include/problem/cec-2017/";

int main() {
  std::unique_ptr<DE::Problem::CECFunction<double>> f;
  std::ofstream file_out("results.csv");
  if (!file_out.good()) {
    std::cout << "Could not open the output CSV file. Exiting." << std::endl;
    return -1;
  }
  file_out << "Function,";
  for (std::size_t i = 0; i < 4; ++i)
    file_out << "Best,Worst,Median,Mean,Std";
  file_out << std::endl;
  file_out.flush();
  for (std::size_t func = 1; func <= 15; ++func) {
    for (const auto& D : {10, 30, 50, 100}) {
      std::vector<double> results;
      for (std::size_t run = 0; run < 52; ++run) {
        initialize_function(f, func, D);
        if (D == 10 && run == 0)
          file_out << f->get_name() << ",";
        DE::Algorithm::SHADE<double> shade(std::move(f), true);  // L-SHADE
        shade.evolve_population(5e10);
        auto solution = shade.get_best();
        results.push_back(solution.best_fitness);
      }
      Stats<double> stats(results);
      auto s = stats.calculate_all();
      file_out << std::scientific << s.min << "," << s.max << "," << s.median
               << "," << std::sqrt(s.variance);
      file_out.flush();
    }
    file_out << std::endl;
  }
  return 0;
}
