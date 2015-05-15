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
  /*
  if (tempabund.size() != tempc.size())
  {
    std::cout << "error in species object when combining two";
  }
  */
  /*
  for (int i = 0 ; i < abundances.size() ; ++i)
  {
    m_traits.push_back(traits[i]);
    m_abundances.push_back(abundances[i]);
  }
  */
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

void Species::Setup(long abund_in, long c_in)
{
  Clear();
  m_abundances.push_back(abund_in);
  m_traits.push_back(c_in);
}
