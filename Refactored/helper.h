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

#ifndef HELPER_H
#define HELPER_H

#include <string>

// a set of tools for converting numbers to letters for species names in phylogenies
std::string alpha_format(long x);
std::string alpha_format_auto(long x_in);
std::string alpha_format(long x_in, int num_digits);

#endif // HELPER_H
