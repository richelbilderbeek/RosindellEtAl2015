#ifndef NDEBUG
#include "ntsim.h"

#include <cassert>

void NTsim::Test() noexcept
{
  {
    static bool is_tested{false};
    if (is_tested) return;
    is_tested = true;
  }
  assert(!"OK");
}

#endif
