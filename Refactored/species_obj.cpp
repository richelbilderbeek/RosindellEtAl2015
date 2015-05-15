#include "species_obj.h"

#include <algorithm>
#include <cassert>

Species::Species()
  : m_traits{},
    m_abundances{}
{

}

void Species::Add(const Species& x) noexcept
{
  std::vector<long> abundances = x.GetAbundances();
  std::vector<long> traits = x.GetTraits();
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
  const double trait_variance = GetTraitVarianceOld();
  const double trait_variance_old = GetTraitVarianceOld();
  assert(trait_variance == trait_variance_old);
  return trait_variance;
}

double Species::GetTraitVarianceOld() const noexcept
{
  if ((m_traits.size() > 1)&&(GetSumAbundances() > 1))
  {
    double tot = 0;
    double tot2 = 0;
    double numr = double(GetSumAbundances());
    for (long i = 0 ; i < m_traits.size() ; ++i)
    {
      tot += double(m_traits[i]*m_abundances[i]);
      tot2 += double(m_traits[i])*double(m_traits[i]*m_abundances[i]);
    }
    return ((numr/(numr-1)) * ((tot2/numr)-((tot/numr)*(tot/numr))));
  }
  else
  {
    return 0;
  }
}

void Species::Setup(long abund_in, long c_in) noexcept
{
  Clear();
  m_abundances.push_back(abund_in);
  m_traits.push_back(c_in);
}
