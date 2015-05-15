#include "species.h"

#include <algorithm>
#include <cassert>
#include <fstream>

Species::Species()
  : m_traits{},
    m_abundances{}
{

}

void Species::Add(const Species& x) noexcept
{
  const auto& abundances = x.GetAbundances();
  const auto& traits = x.GetTraits();
  assert(x.GetAbundances().size() == x.GetTraits().size());
  std::copy(std::begin(abundances),std::end(abundances),std::back_inserter(m_abundances));
  std::copy(std::begin(traits),std::end(traits),std::back_inserter(m_traits));
}

void Species::Clear() noexcept
{
  m_traits.clear();
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

long Species::GetTraitMax() const noexcept
{
  const auto trait_max = *std::max_element(std::begin(m_traits),std::end(m_traits));
  return trait_max;
}

long Species::GetTraitMin() const noexcept
{
  const auto trait_min = *std::min_element(std::begin(m_traits),std::end(m_traits));
  return trait_min;
}

double Species::GetTraitMean() const noexcept
{
  assert(m_traits.size() == m_abundances.size());
  //SUM(trait * abunance)
  const double trait_sum_inner_product =
    std::inner_product(
      std::begin(m_traits),
      std::end(m_traits),
      std::begin(m_abundances),
      0.0
    );
  const double trait_mean
    = trait_sum_inner_product
    / static_cast<double>(GetSumAbundances())
  ;
  return trait_mean;
}

double Species::GetTraitVariance() const noexcept
{
  const double numr = double(GetSumAbundances());

  //TODO: This can be shortened after adding higher-level tests
  double trait_variance = 0.0;
  if (m_traits.size() > 1 && numr > 1)
  {
    //We can calculate a variance
    const double sip = //Sum of Inner Product
      std::inner_product(
        std::begin(m_traits),
        std::end(m_traits),
        std::begin(m_abundances),
        0.0
      );
    const double sipst = //Sum of Inner Product with squared traits
      std::inner_product(
        std::begin(m_traits),
        std::end(m_traits),
        std::begin(m_abundances),
        0.0,
        std::plus<double>(),
        [](const double trait, const double abundance) { return trait * trait * abundance; }
      );

    const double sa = static_cast<double>(GetSumAbundances());
    trait_variance = ((sa/(sa-1.0)) * ((sipst/sa)-((sip/sa)*(sip/sa))));
  }
  return trait_variance;
}

void Species::Setup(long abund_in, long c_in) noexcept
{
  Clear();
  m_abundances.push_back(abund_in);
  m_traits.push_back(c_in);
}
