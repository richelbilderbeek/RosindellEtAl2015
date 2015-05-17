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

#ifndef NDEBUG
#include "ntsim.h"

#include <cassert>
#include <iterator>
#include <iostream>

void NTsim::Test() noexcept
{
  {
    static bool is_tested{false};
    if (is_tested) return;
    is_tested = true;
  }
  const bool verbose{false};
  NTsim sim;
  const int seed = 42;
  const int metacommunity_size = 100;
  const double mutation_rate = 0.1;
  const double selection_strength = 0.1;
  sim.setup(
    seed,
    metacommunity_size,
    mutation_rate,
    selection_strength
  );

  if (verbose)
  {
    const auto& cs = sim.GetFitnessCategories();
    std::copy(std::begin(cs),std::end(cs),
      std::ostream_iterator<long>(std::cout," ")
    );
    std::cout << std::endl;
  }

  for (int i=0; i!=100; ++i)
  {
    for (int j=0; j!=100; ++j)
    {
      sim.sim_step();
    }
    if (verbose)
    {
      const auto& cs = sim.GetFitnessCategories();
      std::copy(std::begin(cs),std::end(cs),std::ostream_iterator<long>(std::cout," "));
      std::cout << std::endl;
    }
  }
  //assert(!"DONE");
}

#endif
