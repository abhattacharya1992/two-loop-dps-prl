// Runs the compiled algebraic part of ggvvamp-1.0 in quad precision on the polylogarithm
// values written by check_cxx.py, and prints the nine computed coefficients A_j at one loop and
// the C_F and C_A parts at two loops, in the q_T scheme.
//
//   ggvv_cxx R_VALUES > A_VALUES
//
// R_VALUES: lines "index value". Entries that are not given are set to NaN, so that a read of
// a value the Python side did not supply shows up as NaN in the output.
#include "ggvvamp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

quadtype harness_CF = 0, harness_CA = 0;

int main(int argc, char** argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s R_VALUES\n", argv[0]);
    return 2;
  }
  ggvvamp<quadtype> g;
  g.r = new quadtype[ggvvamp<quadtype>::rsize];
  for (int i = 0; i < ggvvamp<quadtype>::rsize; ++i) g.r[i] = nanq("");
  FILE* f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 1; }
  char line[512];
  int n = 0;
  while (fgets(line, sizeof line, f)) {
    char* end;
    long i = strtol(line, &end, 10);
    if (end == line) continue;
    if (i <= 0 || i >= ggvvamp<quadtype>::rsize) { fprintf(stderr, "bad index %ld\n", i); return 1; }
    g.r[i] = strtoflt128(end, NULL);
    ++n;
  }
  fclose(f);
  fprintf(stderr, "%d values read\n", n);
#include "calls.inc"
  const int computed[9] = {1, 2, 4, 5, 8, 9, 16, 17, 18};
  const char* part[3] = {"1l", "2l_CF", "2l_CA"};
  for (int p = 1; p <= 2; ++p) {
    harness_CF = (p == 1) ? 1 : 0;
    harness_CA = (p == 2) ? 1 : 0;
    g.p2464();
    for (int k = 0; k < 9; ++k) {
      int j = computed[k];
      char re[64], im[64];
      if (p == 1) {
        quadmath_snprintf(re, sizeof re, "%.34Qe", g.A[j][0][0]);
        quadmath_snprintf(im, sizeof im, "%.34Qe", g.A[j][0][1]);
        printf("%d %s %s %s\n", j, part[0], re, im);
      }
      quadmath_snprintf(re, sizeof re, "%.34Qe", g.A[j][1][0]);
      quadmath_snprintf(im, sizeof im, "%.34Qe", g.A[j][1][1]);
      printf("%d %s %s %s\n", j, part[p], re, im);
    }
  }
  return 0;
}
