#ifndef SPECIES_OBJ_H
#define SPECIES_OBJ_H

#include <iostream>
#include <vector>

///Good species object
struct Species
{
  Species();

  void Add(const Species& x);

  void Clear();

  long GetNumberOfSpecies();

  long GetSumAbundances();
  double GetTraitMean();
  long GetTraitMin();
  long GetTraitMax();
  double GetTraitVariance();

  void Setup(long abund_in, long c_in);

private:

  //E.g. consider a good and incipition species, where
  // both trait and abundances are kept together
  //
  //   |
  //   |
  // +-+-+
  // |   |
  // +   +
  // A   B
  //
  //Each species has a collection of traits
  std::vector<long> m_traits;
  std::vector<long> m_abundances; //Number of individuals

  const std::vector<long>& GetAbundances() const noexcept { return m_abundances; }
  double GetTraitMeanOld();
  long GetTraitMaxOld();
  long GetTraitMinOld();
  double GetTraitVarianceOld();
  const std::vector<long>& GetTraits() const noexcept { return m_traits; }
};

#endif // SPECIES_OBJ_H
