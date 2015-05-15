#ifndef SPECIES_H
#define SPECIES_H

#include <vector>

///Good species object
///E.g. consider a good and incipition species, where
/// both trait and abundances are kept together
/*

   |
   |
 +-+-+
 |   |
 +   +
 A   B

*/
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

  //Each species has a collection of traits
  std::vector<long> m_traits;

  ///Number of individuals
  ///Q: Can one not assume that abundances are above zero?
  ///A: No, abundances can be zero
  std::vector<long> m_abundances;

  const std::vector<long>& GetAbundances() const noexcept { return m_abundances; }
  const std::vector<long>& GetTraits() const noexcept { return m_traits; }
};

#endif // SPECIES_H
