#include "species.h"

#include <algorithm>
#include <cassert>
#include <fstream>

Species::Species()
  : m_fitnesses{},
    m_abundances{}
{

}

void Species::Add(const Species& x) noexcept
{
  const auto& abundances = x.GetAbundances();
  const auto& fitnesses = x.GetFitnesses();
  assert(x.GetAbundances().size() == x.GetFitnesses().size());
  std::copy(std::begin(abundances),std::end(abundances),std::back_inserter(m_abundances));
  std::copy(std::begin(fitnesses),std::end(fitnesses),std::back_inserter(m_fitnesses));
}

void Species::Clear() noexcept
{
  m_fitnesses.clear();
  m_abundances.clear();
}

long Species::GetNumberOfSpecies() const noexcept
{
  //Q: Can one not assume that abundances are above zero?
  //A: No, abundances can be zero
  //for (const auto abundance: m_abundances) { assert(abundance > 0); }

  const auto n = std::count_if(std::begin(m_abundances),std::end(m_abundances),
    [](const auto abundance) { return abundance > 0; }
  );
  return n;
}

long Species::GetSumAbundances() const noexcept
{
  const auto sum = std::accumulate(std::begin(m_abundances),std::end(m_abundances),0);
  return sum;
}

long Species::GetFitnessMax() const noexcept
{
  const auto fitness_max = *std::max_element(std::begin(m_fitnesses),std::end(m_fitnesses));
  return fitness_max;
}

long Species::GetFitnessMin() const noexcept
{
  const auto fitness_min = *std::min_element(std::begin(m_fitnesses),std::end(m_fitnesses));
  return fitness_min;
}

double Species::GetFitnessMean() const noexcept
{
  assert(m_fitnesses.size() == m_abundances.size());
  //SUM(fitness * abunance)
  const double fitness_sum_inner_product =
    std::inner_product(
      std::begin(m_fitnesses),
      std::end(m_fitnesses),
      std::begin(m_abundances),
      0.0
    );
  const double fitness_mean
    = fitness_sum_inner_product
    / static_cast<double>(GetSumAbundances())
  ;
  return fitness_mean;
}

double Species::GetFitnessVariance() const noexcept
{
  const double numr = double(GetSumAbundances());

  //TODO: This can be shortened after adding higher-level tests
  double fitness_variance = 0.0;
  if (m_fitnesses.size() > 1 && numr > 1)
  {
    //We can calculate a variance
    const double sip = //Sum of Inner Product
      std::inner_product(
        std::begin(m_fitnesses),
        std::end(m_fitnesses),
        std::begin(m_abundances),
        0.0
      );
    const double sipst = //Sum of Inner Product with squared fitnesses
      std::inner_product(
        std::begin(m_fitnesses),
        std::end(m_fitnesses),
        std::begin(m_abundances),
        0.0,
        std::plus<double>(),
        [](const double fitness, const double abundance) { return fitness * fitness * abundance; }
      );

    const double sa = static_cast<double>(GetSumAbundances());
    fitness_variance = ((sa/(sa-1.0)) * ((sipst/sa)-((sip/sa)*(sip/sa))));
  }
  return fitness_variance;
}

void Species::Setup(long abund_in, long c_in) noexcept
{
  Clear();
  m_abundances.push_back(abund_in);
  m_fitnesses.push_back(c_in);
}
