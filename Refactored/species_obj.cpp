#include "species_obj.h"

#include <algorithm>
#include <cassert>

Species::Species()
  : m_traits{},
    m_abundances{}
{

}

void Species::Add(const Species& x)
{
  std::vector<long> abundances = x.GetAbundances();
  std::vector<long> traits = x.GetTraits();
  assert(x.GetAbundances().size() == x.GetTraits().size());
  std::copy(std::begin(abundances),std::end(abundances),std::back_inserter(m_abundances));
  std::copy(std::begin(traits),std::end(traits),std::back_inserter(m_traits));
}

void Species::Clear()
{
  m_traits.clear();
  m_abundances.clear();
}


long Species::GetNumberOfSpecies()
{
  //TODO: Can one not assume that abundances are above zero?

  const auto n = std::count_if(
    std::begin(m_abundances),
    std::end(m_abundances),
    [](const auto abundance) { return abundance > 0; }
  );
  return n;
}

long Species::GetSumAbundances()
{
  const auto sum = std::accumulate(std::begin(m_abundances),std::end(m_abundances),0);
  return sum;
}

long Species::GetTraitMax()
{
  const auto trait_max = GetTraitMaxOld();
  const auto trait_max_old = GetTraitMaxOld();
  assert(trait_max == trait_max_old);
  return trait_max;
}

long Species::GetTraitMaxOld()
{
  if (m_traits.size() > 0)
  {
    long toret = m_traits[0];
    for (long i = 1 ; i < m_traits.size() ; ++i)
    {
      if (toret < m_traits[i])
      {
        toret = m_traits[i];
      }
    }
    return toret;
  }
  else
  {
    return 0;
  }
}

long Species::GetTraitMin()
{
  const auto trait_min = GetTraitMinOld();
  const auto trait_min_old = GetTraitMinOld();
  assert(trait_min == trait_min_old);
  return trait_min;
}

long Species::GetTraitMinOld()
{
  if (m_traits.size() > 0)
  {
    long toret = m_traits[0];
    for (long i = 1 ; i < m_traits.size() ; ++i)
    {
      if (toret > m_traits[i])
      {
        toret = m_traits[i];
      }
    }
    return toret;
  }
  else
  {
    return 0;
  }
}

double Species::GetTraitMean()
{
  const double trait_mean = GetTraitMeanOld();
  const double trait_mean_old = GetTraitMeanOld();
  assert(trait_mean == trait_mean_old);
  return trait_mean;
}

double Species::GetTraitMeanOld()
{
  if ((m_traits.size() > 0)&&(GetSumAbundances() > 0))
  {
    double toret = 0;
    for (long i = 0 ; i < m_traits.size() ; ++i)
    {
      toret += double(m_traits[i]*m_abundances[i]);
    }
    return (toret/double(GetSumAbundances()));
  }
  else
  {
    return 0;
  }
}

double Species::GetTraitVariance()
{
  const double trait_variance = GetTraitVarianceOld();
  const double trait_variance_old = GetTraitVarianceOld();
  assert(trait_variance == trait_variance_old);
  return trait_variance;
}

double Species::GetTraitVarianceOld()
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

void Species::Setup(long abund_in, long c_in)
{
  Clear();
  m_abundances.push_back(abund_in);
  m_traits.push_back(c_in);
}
