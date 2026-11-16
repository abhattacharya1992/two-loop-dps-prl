// Compare the original GiNaC special-function stage without Python inputs.
// Build with the unmodified ggvvamp.cpp and p0002.cpp from ggvvamp-1.0.
#include <ggvvamp.h>
#include <ginac/ginac.h>
#include <cln/float.h>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  if (argc != 6) throw std::runtime_error("usage: native_functions DIGITS S T MA2 MB2");
  using GiNaC::numeric;
  using GiNaC::ex_to;
  GiNaC::Digits = std::stoi(argv[1]);
  cln::default_float_format = cln::float_format(std::stoi(argv[1]));
  const numeric s(argv[2]), t(argv[3]), ma2(argv[4]), mb2(argv[5]);
  const numeric u = ma2 + mb2 - s - t;
  const numeric root = sqrt((s-ma2-mb2)*(s-ma2-mb2)-4*ma2*mb2);
  const numeric t0=(ma2+mb2-s-root)/2, t1=(ma2+mb2-s+root)/2;
  if (ma2<=0 || mb2<=0 || t<t0 || t>t1) throw std::runtime_error("nonphysical point");
  const numeric vars[] = {ma2,t/ma2,u/ma2,mb2/ma2,-t0/ma2,
                         t1*(u+t)/(ma2*mb2)-1,t1*t/(ma2*mb2)};
  ggvvamp<cln::cl_F> g(5,IRschemeQt);
  for (int k=0;k<7;++k) g.r[k+1]=ggvvamp<cln::cl_F>::n(vars[k]);
  g.p0002();
  for (int k=1;k<=1620;++k) std::cout << "R " << k << " " << g.r[k] << "\n";
}
