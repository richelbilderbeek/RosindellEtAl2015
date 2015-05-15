#ifndef HELPER_H
#define HELPER_H

#include <string>

// a set of tools for converting numbers to letters for species names in phylogenies
std::string alpha_format(long x);
std::string alpha_format_auto(long x_in);
std::string alpha_format(long x_in, int num_digits);

#endif // HELPER_H
