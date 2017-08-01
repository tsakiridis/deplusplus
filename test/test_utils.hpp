#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <vector>
#include <string>

extern std::vector<std::vector<double>> x_tests;

void init_values();

std::string shift_file(const std::size_t problem);
std::string rotation_file(const std::size_t problem, const std::size_t D);
std::string shuffle_file(const std::size_t problem, const std::size_t D);
std::string bias_file(const std::size_t problem);

void load_shift_data(const int nx, const int func_num);
void load_shuffle_data(const int nx, const int cf_num, const int func_num);
void load_bias_data(const int cf_num, const int func_num);

#endif
