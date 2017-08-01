#include "test_utils.hpp"
#include <ctime>
#include <vector>
#include <algorithm>
#include <string>
#include "rand.hpp"
#include "cec17_test_func.hpp"

// rand.hpp
size_t SEED = std::time(nullptr);

// Test values
constexpr std::size_t number_of_tests = 25;
const std::vector<std::size_t> allowed_ranges = {10, 30, 50, 100};
std::vector<std::vector<double>> x_tests(number_of_tests);

void init_values() {
  for (auto& vec : x_tests) {
    vec.resize(allowed_ranges[rand_uniform_int(0, allowed_ranges.size() - 1)]);
    for (auto& val : vec)
      val = rand_uniform_real(-100.0, 100.0);
  }
  const std::size_t max_size =
      std::max_element(x_tests.begin(), x_tests.end(), [](const auto& vec_1,
                                                          const auto& vec_2) {
        return vec_1.size() < vec_2.size();
      })->size();
  y = new double[max_size];
  z = new double[max_size];
}

// Shift and rotation files
const std::string base_path = "cec-2017/";

std::string shift_file(const std::size_t problem) {
  return base_path + "shift_data_" + std::to_string(problem) + ".txt";
}

std::string rotation_file(const std::size_t problem, const std::size_t D) {
  return base_path + "M_" + std::to_string(problem) + "_D" + std::to_string(D) +
         ".txt";
}

std::string shuffle_file(const std::size_t problem, const std::size_t D) {
  return base_path + "shuffle_data_" + std::to_string(problem) + "_D" +
         std::to_string(D) + ".txt";
}

void load_shift_data(const int nx, const int func_num) {
  char FileName[256], tmpchar;
  const int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
  int cf_num = cf_nums[func_num];
  /* Load shift_data */
  std::string format = base_path + "shift_data_%d.txt";
  sprintf(FileName, format.c_str(), func_num);
  FILE* fpt = fopen(FileName, "r");
  if (fpt == NULL) {
    printf("\n Error: Cannot open input file %s for reading \n", FileName);
  }
  OShift = (double*)malloc(cf_num * nx * sizeof(double));
  if (OShift == NULL) {
    printf("\nError: there is insufficient memory available!\n");
    exit(1);
  }

  for (int i = 0; i < nx * cf_nums[func_num]; i++) {
    fscanf(fpt, "%lf", &OShift[i]);
    // if finish reading a row, skip the rest of the line to the next line
    // bug fixed 2014/12/26 by Q. Chen
    if (cf_nums[func_num] > 1 && ((i + 1) % nx) == 0) {
      // printf("%d\t%d\t%d\n", cf_nums[func_num], nx, i);
      fscanf(fpt, "%c", &tmpchar);
      while (tmpchar != '\n') {
        // printf("%c", tmpchar);
        fscanf(fpt, "%c", &tmpchar);
      }
    }
  }
  fclose(fpt);
}

void load_shuffle_data(const int nx, const int cf_num, const int func_num) {
  const int bShuffle[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0};
  /* Load Shuffle_data */
  if (bShuffle[func_num] == 1) {
    char FileName[256];
    std::string format = base_path + "shuffle_data_%d_D%d.txt";
    sprintf(FileName, format.c_str(), func_num, nx);
    FILE* fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      printf("\n Error: Cannot open input file for reading \n");
      exit(1);
    }
    SS = (int*)malloc(nx * cf_num * sizeof(int));
    if (SS == NULL) {
      printf("\nError: there is insufficient memory available!\n");
      exit(1);
    }
    for (int i = 0; i < nx * cf_num; i++) {
      fscanf(fpt, "%d", &SS[i]);
    }
    fclose(fpt);
  }
}
