#ifndef SPECIES_OBJ_H
#define SPECIES_OBJ_H

#include <iostream>
#include <vector>

///Good species object
struct Species
{
  Species();

  void Add(const Species& x) noexcept;

  void Clear() noexcept;

  long GetNumberOfSpecies() const noexcept;
  long GetSumAbundances() const noexcept;
  double GetTraitMean() const noexcept;
  long GetTraitMin() const noexcept;
  long GetTraitMax() const noexcept;
  double GetTraitVariance() const noexcept;

  void Setup(long abund_in, long c_in) noexcept;

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
  double GetTraitMeanOld() const noexcept;
  long GetTraitMaxOld() const noexcept;
  long GetTraitMinOld() const noexcept;
  double GetTraitVarianceOld() const noexcept;
  const std::vector<long>& GetTraits() const noexcept { return m_traits; }
};

#endif // SPECIES_OBJ_H
