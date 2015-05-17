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

#include "nrrand.h"

NRrand::NRrand()
  : m_idum{0},
    m_j{0},
    m_k{0},
    m_idum2{0},
    m_iy{0},
    m_iv{},
    m_seeded{false},
    m_temp{0.0}
{

}


