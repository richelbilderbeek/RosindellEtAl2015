// Code for the article
//
//   Rosindell, James, Luke J. Harmon, and Rampal S. Etienne.
//   "Unifying ecology and macroevolution with individual‐based theory."
//   Ecology letters 18.5 (2015): 472-482.
//
// Original version:
//   Version 1.2 Jan 2014
//   Author: James Rosindell
//
// Refactoring by Richel Bilderbeek

#ifndef SPECIES_H
#define SPECIES_H

#include <vector>

///Good species object
///E.g. consider a good and incipition species, where
/// both fitness category (symbol used in article: c)
/// and abundances are kept together
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
  double GetFitnessMean() const noexcept;
  long GetFitnessMin() const noexcept;
  long GetFitnessMax() const noexcept;
  double GetFitnessVariance() const noexcept;
  void Setup(const long abundance, const long fitness) noexcept;

private:

  ///Each species has a collection of fitness categories
  ///Symbol used for fitness category: c
  std::vector<long> m_fitnesses;

  ///Number of individuals
  ///Q: Can one not assume that abundances are above zero?
  ///A: No, abundances can be zero
  std::vector<long> m_abundances;

  const std::vector<long>& GetAbundances() const noexcept { return m_abundances; }
  const std::vector<long>& GetFitnesses() const noexcept { return m_fitnesses; }
};

#endif // SPECIES_H
