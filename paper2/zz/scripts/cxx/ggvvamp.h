// Minimal replacement for ggvvamp.h of ggvvamp-1.0, used only by check_cxx.py.
//
// It lets the generated C++ files p0003.cpp ... p2464.cpp of von Manteuffel and Tancredi be
// compiled unchanged (apart from the Casimirs in p2464.cpp, see check_cxx.py) with the quad
// precision type of the original code, without GiNaC or CLN. The polylogarithms of p0002.cpp
// are not computed here: their values are read from a file written by the Python translation,
// so that the check covers the algebraic part (p0003-p2463) and the assignment (p2464).
#ifndef GGVVAMP_H
#define GGVVAMP_H

#include <stdint.h>

#define QUAD_PRECISION
extern "C" {
#include <quadmath.h>
}
typedef __float128 quadtype;

namespace GiNaC {
  struct pi_tag {};
  static const pi_tag Pi = pi_tag();
}

const int nAcoeff = 20; // coefficient type {1,...,20}
const int nloop = 2;    // {0=1-loop, 1=2-loop}
const int nreim = 2;    // {0=Re(), 1=Im()}

typedef int64_t longint;

extern quadtype harness_CF, harness_CA;   // the Casimirs used by the patched p2464.cpp

template<class numtype>
class ggvvamp {
public:
  numtype A[nAcoeff+1][nloop][nreim];
  numtype* r;
  // as in the original header for quadtype
  static inline numtype n(longint a) { return static_cast<numtype>(a); }
  static inline numtype n(longint a, longint b) { return a / static_cast<numtype>(b); }
  static inline numtype n(GiNaC::pi_tag) { return M_PIq; }
#include "ggvvamp.xh"
};

#endif
