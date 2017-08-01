/*!
 * This example demonstrates the distributed use of DE++ (using a thread pool)
 * for CEC benchmark functions.
 *
 * If you haven't done so already, you need to generate the shift, rotation
 * etc. files using the Matlab / Octave script. Kindly check README.md.
 *
 * To define the base path where these files are located you need to specify
 * it using the "-p" flag.
 *
 * Additionally you can select a new SEED using the "-s" flag.
 *
 * This example uses L-SHADE and runs all experiments; it further accumulates
 * the results and writes them out to a csv file.
 */

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <unistd.h>

#include "problem/cec_all_functions.hpp"
#include "algorithm/shade.hpp"
#include "stats.hpp"
#include "ThreadPool.h"

size_t SEED = 0;
std::string base_path;

inline bool ends_with(const std::string& value, const std::string& ending) {
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

void write_statistics(std::ofstream& log, const std::vector<double>& samples) {
  Stats<double> stats(samples);
  auto s = stats.calculate_all();
  log << std::scientific << "," << s.min << "," << s.max << "," << s.median
      << "," << s.mean << "," << std::sqrt(s.variance);
}

struct ResultsOneRun {
  const char* name;
  const std::size_t function_id;
  const std::size_t dimensions;
  const std::size_t run;
  const double fitness;

  ResultsOneRun(const char* n,
                const std::size_t f_id,
                const std::size_t d,
                const std::size_t r,
                const double f)
      : name(n), function_id(f_id), dimensions(d), run(r), fitness(f) {}
};

int main(int argc, char** argv) {
  // Parse command-line arguments
  int c;
  while ((c = getopt(argc, argv, "p:s:")) != -1)
    switch (c) {
      case 'p':
        base_path = std::string(optarg);
        break;
      case 's':
        SEED = size_t(optarg);
        break;
      case '?':
        std::cout << "Unkown parameter: " << optopt << std::endl;
        return -1;
    }
  if (base_path.empty()) {
    std::cout << "You need to specify a base path!" << std::endl;
    return -1;
  } else if (!ends_with(base_path, "/")) {
    base_path.push_back('/');
  }

  // Generate the thread pool
  ThreadPool pool(std::thread::hardware_concurrency());
  std::vector<std::future<ResultsOneRun>> results;
  for (std::size_t func = 14; func <= 14; ++func) {
    for (const auto& D : {10}) {  //, 30, 50, 100}) {
      for (std::size_t run = 0; run < 52; ++run) {
        results.emplace_back(pool.enqueue([func, D, run] {
          std::unique_ptr<DE::Problem::CECFunction<double>> f;
          initialize_function(f, func, D);
          auto name = f->get_name();
          DE::Algorithm::SHADE<double> shade(std::move(f), true);  // L-SHADE
          shade.evolve_population(5e8);
          auto solution = shade.get_best();
          return ResultsOneRun(name, func, D, run, solution.best_fitness);
        }));
      }
    }
  }

  // Accumulate and log results
  std::ofstream file_out("results10.csv");
  file_out << "Function";
  for (std::size_t i = 0; i < 4; ++i)
    file_out << ",Best,Worst,Median,Mean,Std";
  std::size_t function_num = 0, current_dim = 10;
  std::vector<double> samples;
  for (auto&& result : results) {
    const auto& res = result.get();
    if (current_dim != res.dimensions) {
      current_dim = res.dimensions;
      write_statistics(file_out, samples);
      samples.clear();
    }
    samples.push_back(res.fitness);
    if (res.function_id > function_num) {
      file_out << std::endl << res.name;
      ++function_num;
    }
  }
  write_statistics(file_out, samples);  // final results
  return 0;
}
