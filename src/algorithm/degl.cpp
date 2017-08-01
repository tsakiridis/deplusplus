#include "algorithm/degl.hpp"
#include "algorithm/differential_evolution.hpp"
#include "rand.hpp"

namespace DE {
namespace Algorithm {

template <class T>
constexpr float DEGL<T>::F;

template <class T>
constexpr float DEGL<T>::Cr;

template <class T>
DEGL<T>::DEGL(std::shared_ptr<Problem::Base<T>> problem,
              const std::vector<T>& initial_chromosome,
              const bool minimize)
    : Base<T>(problem,
              initial_chromosome,
              10 * problem->get_number_of_genes(),
              minimize),
      k_(Base<T>::N_ <= 10 ? 1 : (Base<T>::N_ / 10) + 1) {
  initialize_weights();
}

template <class T>
DEGL<T>::DEGL(std::shared_ptr<Problem::Base<T>> problem, const bool minimize)
    : Base<T>(problem, 10 * problem->get_number_of_genes(), minimize),
      k_(Base<T>::N_ <= 10 ? 1 : (Base<T>::N_ / 10) + 1) {
  initialize_weights();
}

template <class T>
void DEGL<T>::evolve_population(const std::size_t max_generations) {
#ifdef CEC_MAX_EVALUATIONS
  const std::size_t max_evaluations = Base<T>::D_ * 10e4;
  std::size_t evaluations = 0;
#endif
  for (std::size_t g = 0; g < max_generations; ++g) {
    std::size_t best_index = Base<T>::best_index();
    for (std::size_t i = 0; i < Base<T>::N_; ++i) {
      auto u = mutate(i, best_index);
      auto v = binary_crossover(Base<T>::x_[i], u, Cr);
      auto trial_fitness = Base<T>::p_problem_->fitness(v);
      if (Base<T>::compare_fitnesses_with_equality(Base<T>::fit_[i],
                                                   trial_fitness)) {
        std::swap(Base<T>::x_[i], v);
        Base<T>::fit_[i] = trial_fitness;
        w_[i] = w_mutated_[i];
      }
    }
#ifdef CEC_MAX_EVALUATIONS
    evaluations += Base<T>::N_;
    if (evaluations > max_evaluations)
      break;
#endif
  }
}

/*! Private member functions */

template <class T>
void DEGL<T>::initialize_weights() {
  w_.resize(Base<T>::N_);
  w_mutated_.resize(Base<T>::N_);
  for (auto& w : w_)
    w = rand_uniform_real(0.05, 0.95);
}

template <class T>
std::vector<std::size_t> DEGL<T>::get_neighborhood_indexes(
    const std::size_t index) const {
  std::vector<std::size_t> places;
  places.reserve(2 * k_);
  // Forward
  std::size_t fwd = index, counter = 0;

  while (++counter <= k_) {
    if (++fwd == Base<T>::N_)
      fwd = 0;

    places.push_back(fwd);
  }

  // Backwards
  counter = 0;
  std::size_t bck = index;

  while (++counter <= k_) {
    if (bck-- == 0)
      bck = Base<T>::N_ - 1;

    places.push_back(bck);
  }

  return places;
}

template <class T>
std::array<std::size_t, 2> DEGL<T>::get_random_local_indexes(
    const std::vector<std::size_t>& places) const {
  std::array<std::size_t, 2> rand;
  std::size_t r1, r2;
  r1 = rand_uniform_int(0, places.size() - 1);
  rand[0] = places[r1];

  do {
    r2 = rand_uniform_int(0, places.size() - 1);
  } while (r2 == r1);

  rand[1] = places[r2];
  return rand;
}

template <class T>
std::array<std::size_t, 2> DEGL<T>::get_random_global_indexes(
    const std::size_t avoid) const {
  std::array<std::size_t, 2> r;
  do {
    r[0] = rand_uniform_int(0, Base<T>::N_ - 1);
  } while (r[0] == avoid);
  do {
    r[1] = rand_uniform_int(0, Base<T>::N_ - 1);
  } while (r[1] == r[0] || r[1] == avoid);
  return r;
}

template <class T>
std::size_t DEGL<T>::find_local_best(
    const std::vector<std::size_t>& places) const {
  double best = Base<T>::fit_[places[0]];
  std::size_t pbest = places[0];

  for (std::size_t i = 1; i < places.size(); ++i) {
    if (Base<T>::compare_fitnesses(best, Base<T>::fit_[places[i]])) {
      best = Base<T>::fit_[places[i]];
      pbest = places[i];
    }
  }

  return pbest;
}

template <class T>
std::vector<T> DEGL<T>::mutate(const std::size_t index,
                               const std::size_t global_best) {
  const auto places = get_neighborhood_indexes(index);
  const auto local_rand = get_random_local_indexes(places);
  const auto local_best = find_local_best(places);
  const auto global_rand = get_random_global_indexes(index);

  std::vector<T> v(Base<T>::x_[index]);
  w_mutated_[index] = w_[index] + F * (w_[global_best] - w_[index]) +
                      F * (w_[global_rand[0]] - w_[global_rand[1]]);
  w_mutated_[index] = std::max(0.05, std::min(0.95, double(w_mutated_[index])));
  for (std::size_t i = 0; i < Base<T>::D_; ++i) {
    v[i] =
        w_mutated_[index] *
            (this->x_[index][i] +
             F * (this->x_[global_best][i] - this->x_[index][i]) +
             F * (this->x_[global_rand[0]][i] - this->x_[global_rand[1]][i])) +
        (1 - w_mutated_[index]) *
            (this->x_[index][i] +
             F * (this->x_[local_best][i] - this->x_[index][i]) +
             F * (this->x_[local_rand[0]][i] - this->x_[local_rand[1]][i]));
  }
  Base<T>::p_problem_->constrain(v);
  return v;
}

// explicit instantiations
template class DEGL<float>;
template class DEGL<double>;

}  // namespace Algorithm
}  // namespace DE
