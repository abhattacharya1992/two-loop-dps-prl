// GiNaC 1.8 uses append/initializer lists instead of the old lst(a,b) API.
#include <ginac/ginac.h>
template<class A, class B>
GiNaC::lst legacy_lst(const A& a,const B& b) {
  GiNaC::lst result;
  result.append(a);
  result.append(b);
  return result;
}
