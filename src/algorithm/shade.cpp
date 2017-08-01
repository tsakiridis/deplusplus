#include "algorithm/shade.hpp"
#include <numeric>
#include <assert.h>
#include <thread>
#include "algorithm/differential_evolution.hpp"
#include "rand.hpp"

namespace DE {
namespace Algorithm {

template <class T>
SHADE<T>::SHADE(std::shared_ptr<Problem::Base<T>> problem,
                const std::vector<T>& initial_chromosome,
                const bool use_linear_size_reduction,
                const bool minimize)
    : Base<T>(problem,
              initial_chromosome,
              18 * problem->get_number_of_genes(),
              minimize),
      N_(Base<T>::N_),
      p_(0.11 * N_),
      H_(6),
      A_size_(2.6 * N_),
      use_linear_size_reduction_(use_linear_size_reduction) {
  A_.reserve(A_size_);
  top_p_.resize(p_);
  Cr_.resize(H_, 0.5);
  F_.resize(H_, 0.5);
}

template <class T>
SHADE<T>::SHADE(std::shared_ptr<Problem::Base<T>> problem,
                const std::vector<T>& initial_chromosome,
                const double initial_fitness,
                const bool use_linear_size_reduction,
                const bool minimize)
    : Base<T>(problem,
              initial_chromosome,
              18 * problem->get_number_of_genes(),
              minimize),
      N_(Base<T>::N_),
      p_(0.11 * N_),
      H_(6),
      A_size_(2.6 * N_),
      use_linear_size_reduction_(use_linear_size_reduction) {
  A_.reserve(A_size_);
  top_p_.resize(p_);
  Cr_.resize(H_, 0.5);
  F_.resize(H_, 0.5);
  Base<T>::fit_[0] = initial_fitness;
}

template <class T>
SHADE<T>::SHADE(std::shared_ptr<Problem::Base<T>> problem,
                const bool use_linear_size_reduction,
                const bool minimize)
    : Base<T>(problem,
              18 * problem->get_number_of_genes(),
              minimize),
      N_(Base<T>::N_),
      p_(0.11 * N_),
      H_(6),
      A_size_(2.6 * N_),
      use_linear_size_reduction_(use_linear_size_reduction) {
  A_.reserve(A_size_);
  top_p_.resize(p_);
  Cr_.resize(H_, 0.5);
  F_.resize(H_, 0.5);
}

template <class T>
void SHADE<T>::evolve_population(const std::size_t max_generations) {
#ifdef CEC_MAX_EVALUATIONS
  const std::size_t max_evaluations = Base<T>::D_ * 10e4;
  std::size_t evaluations = 0;
#endif
  std::vector<float> S_Cr, S_F;
  std::vector<double> delta_fit;
  const size_t num_threads = std::thread::hardware_concurrency();
  std::vector<std::thread> threads(num_threads);
  for (std::size_t g = 0; g < max_generations; ++g) {
    S_Cr.clear();
    S_F.clear();
    delta_fit.clear();
    update_top_p_solutions();

    if (Base<T>::allow_parallel_) {
      for (size_t i = 0; i < num_threads; ++i)
        threads[i] =
            std::thread(&SHADE<T>::parallel_evolve, this, i, std::ref(S_Cr),
                        std::ref(S_F), std::ref(delta_fit));

      for (size_t i = 0; i < num_threads; ++i)
        threads[i].join();
    } else {
      for (std::size_t i = 0; i < N_; ++i) {
        const auto r_i = rand_uniform_int(0, H_ - 1);
        const auto Cr = get_crossover_factor(r_i);
        const auto F = get_scale_factor(r_i);
        const auto mutant = mutate(i, F);
        auto trial = binary_crossover<T>(Base<T>::x_[i], mutant, Cr);
        const auto trial_fitness = Base<T>::p_problem_->fitness(trial);
        if (Base<T>::compare_fitnesses_with_equality(Base<T>::fit_[i],
                                                     trial_fitness)) {
          std::swap(Base<T>::x_[i], trial);
          if (Base<T>::fit_[i] != trial_fitness) {
            add_to_archive(trial);
            S_Cr.push_back(Cr);
            S_F.push_back(F);
            delta_fit.push_back(fabs(Base<T>::fit_[i] - trial_fitness));
            Base<T>::fit_[i] = trial_fitness;
          }
        }
      }
    }

#ifdef CEC_MAX_EVALUATIONS
    evaluations += N_;
    if (evaluations > max_evaluations)
      break;
#endif
    memory_update(S_Cr, S_F, delta_fit, g % H_);
    if (use_linear_size_reduction_)
      linear_size_reduction(g, max_generations);
  }
}

template <class T>
void SHADE<T>::parallel_evolve(const std::size_t thread_id,
                               std::vector<float>& S_Cr,
                               std::vector<float>& S_F,
                               std::vector<double>& delta_fit) {
  for (std::size_t i = thread_id; i < N_;
       i += std::thread::hardware_concurrency()) {
    const auto r_i = rand_uniform_int(0, H_ - 1);
    pop_mutex_.lock();
    const auto Cr = get_crossover_factor(r_i);
    const auto F = get_scale_factor(r_i);
    const auto mutant = mutate(i, F);
    auto trial = binary_crossover<T>(Base<T>::x_[i], mutant, Cr);
    pop_mutex_.unlock();
    const auto trial_fitness = Base<T>::p_problem_->fitness(trial);
    if (Base<T>::compare_fitnesses_with_equality(Base<T>::fit_[i],
                                                 trial_fitness)) {
      pop_mutex_.lock();
      std::swap(Base<T>::x_[i], trial);
      if (Base<T>::fit_[i] != trial_fitness) {
        add_to_archive(trial);
        S_Cr.push_back(Cr);
        S_F.push_back(F);
        delta_fit.push_back(fabs(Base<T>::fit_[i] - trial_fitness));
        Base<T>::fit_[i] = trial_fitness;
      }
      pop_mutex_.unlock();
    }
  }
}

template <class T>
float SHADE<T>::get_crossover_factor(const std::size_t rand_index) const {
  auto Cr = (Cr_[rand_index] == TERMINAL_VALUE) ? 0 : rand_normal(
                                                          Cr_[rand_index], 0.1);
  return std::min(std::max(Cr, 0.0), 1.0);
}

template <class T>
float SHADE<T>::get_scale_factor(const std::size_t rand_index) const {
  float F;
  do {
    F = rand_cauchy(F_[rand_index], 0.1);
    if (F > 1)
      return 1.0;
  } while (F <= 0.0);
  return F;
}

template <class T>
void SHADE<T>::memory_update(const std::vector<float>& S_Cr,
                             const std::vector<float>& S_F,
                             const std::vector<double>& delta_fit,
                             const std::size_t k) {
  assert(S_Cr.size() == S_F.size());
  assert(S_Cr.size() == delta_fit.size());
  assert(k < H_);
  if (!S_Cr.empty()) {
    std::vector<double> weights(delta_fit);
    auto sum = std::accumulate(delta_fit.begin(), delta_fit.end(), 0.0);
    for (auto& w : weights)
      w /= sum;
    if (Cr_[k] == TERMINAL_VALUE ||
        *std::max_element(S_Cr.begin(), S_Cr.end()) == 0) {
      Cr_[k] = TERMINAL_VALUE;
    } else {
      Cr_[k] = weighted_lehmer_mean(weights, S_Cr);
    }
    F_[k] = weighted_lehmer_mean(weights, S_F);
  }
}

template <class T>
float SHADE<T>::weighted_lehmer_mean(const std::vector<double>& weights,
                                     const std::vector<float>& values) const {
  float mean_num = 0.0, mean_den = 0.0;
  for (std::size_t i = 0; i < weights.size(); ++i) {
    mean_num += weights[i] * values[i] * values[i];
    mean_den += weights[i] * values[i];
  }
  return mean_num / mean_den;
}

template <class T>
void SHADE<T>::add_to_archive(const std::vector<T>& chromosome) {
  if (A_.size() >= A_size_) {  // overwrite random
    std::size_t rand_index = rand_uniform_int(0, A_.size() - 1);
    A_[rand_index] = chromosome;
  } else {
    A_.push_back(chromosome);
  }
}

template <class T>
void SHADE<T>::update_top_p_solutions() {
  std::vector<std::size_t> indices(N_, 0);
  std::iota(indices.begin(), indices.end(), 0);
  std::partial_sort(indices.begin(), indices.begin() + p_, indices.end(),
                    [this](const auto& first, const auto& second) {
                      return !this->compare_fitnesses(this->fit_[first],
                                                      this->fit_[second]);
                    });
  for (std::size_t i = 0; i < p_; ++i)
    top_p_[i] = indices[i];
}

template <class T>
std::vector<T> SHADE<T>::mutate(const std::size_t base_index,
                                const float F) const {
  assert(base_index < N_);
  assert(F > 0);
  auto mutant = Base<T>::x_[base_index];
  std::size_t rand_pbest_index =
                  (p_ <= 1) ? top_p_[0] : top_p_[rand_uniform_int(0, p_ - 1)],
              rand_1, rand_2 = rand_uniform_int(0, N_ + A_.size() - 2);
  do {
    rand_1 = rand_uniform_int(0, N_ - 1);
  } while (rand_1 == base_index);
  const std::vector<T>& x_r_2 =
      (rand_2 >= N_) ? A_[rand_2 - N_] : Base<T>::x_[rand_2];
  for (std::size_t j = 0; j < Base<T>::D_; ++j) {
    mutant[j] += F * (Base<T>::x_[rand_pbest_index][j] - mutant[j]) +
                 F * (Base<T>::x_[rand_1][j] - x_r_2[j]);
  }
  Base<T>::p_problem_->constrain(mutant);
  return mutant;
}

template <class T>
void SHADE<T>::linear_size_reduction(const std::size_t current_generation,
                                     const std::size_t max_generations) {
  std::size_t N =
      (4.0 - 18 * Base<T>::D_) * current_generation / max_generations +
      18 * Base<T>::D_ + 1;
  if (N < N_) {
    std::vector<std::size_t> indices(N_, 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [this](const auto& first,
                                                     const auto& second) {
      return this->compare_fitnesses(this->fit_[second], this->fit_[first]);
    });
    auto new_population = Base<T>::x_;
    std::transform(indices.begin(), indices.end(), new_population.begin(),
                   [this](const auto& i) { return this->x_[i]; });
    std::swap(Base<T>::x_, new_population);
    auto new_fitnesses = Base<T>::fit_;
    std::transform(indices.begin(), indices.end(), new_fitnesses.begin(),
                   [this](const auto& i) { return this->fit_[i]; });
    std::swap(Base<T>::fit_, new_fitnesses);
    N_ = N;
    Base<T>::x_.resize(N_);
    Base<T>::fit_.resize(N_);
    A_size_ = 2.6 * N_;
    if (A_.size() > A_size_)
      A_.resize(A_size_);
    p_ = std::max(1.0, 0.11 * N_);
    top_p_.resize(p_);
  }
}

// explicit instantiations
template class SHADE<float>;
template class SHADE<double>;

}  // namespace Algorithm
}  // namespace DE
