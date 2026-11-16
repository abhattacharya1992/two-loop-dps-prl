// Independent native reader for FORM expressions and bare dimensional residues.
#include <ginac/ginac.h>
#include <cln/cln.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
using C=cln::cl_N;
using Key=std::pair<int,int>;
using Poly=std::map<Key,C>;
C one(){return cln::cl_float(1,cln::default_float_format);}
Poly scalar(C v){return {{{0,0},v}};}
Poly add(Poly a,const Poly&b,int sign=1){for(auto&v:b)a[v.first]=a.count(v.first)?a[v.first]+sign*v.second:sign*v.second;return a;}
Poly mul(const Poly&a,const Poly&b){Poly r;for(auto&x:a)for(auto&y:b){Key k{x.first.first+y.first.first,x.first.second+y.first.second};C v=x.second*y.second;r[k]=r.count(k)?r[k]+v:v;}return r;}
C value(const Poly&a){if(a.size()!=1||a.begin()->first!=Key{0,0})throw std::runtime_error("non-scalar argument");return a.begin()->second;}
Poly divp(Poly a,const Poly&b){if(b.size()!=1)throw std::runtime_error("non-monomial denominator");auto v=*b.begin();Poly r;for(auto&x:a)r[{x.first.first-v.first.first,x.first.second-v.first.second}]=x.second/v.second;return r;}
Poly power(Poly a,const Poly&b){C c=value(b);if(!cln::zerop(cln::imagpart(c)))throw std::runtime_error("complex exponent");auto n=cln::truncate1(cln::realpart(c));if(c!=C(n))throw std::runtime_error("nonintegral exponent");int k=cln::cl_I_to_int(n);if(k<0)return divp(scalar(one()),power(a,scalar(C(-n))));Poly r=scalar(one());while(k){if(k&1)r=mul(r,a);k>>=1;if(k)a=mul(a,a);}return r;}
struct Parser {
 const std::string&s;size_t i=0;std::map<std::string,Poly> symbols;std::map<std::string,C> functions,numeric_functions;
 Parser(const std::string&t,const std::map<std::string,Poly>&v):s(t),symbols(v){}
 void ws(){while(i<s.size()&&isspace(s[i]))++i;}
 [[noreturn]]void fail(const std::string&m){throw std::runtime_error(m+" at "+std::to_string(i)+" near "+s.substr(i,60));}
 int prec(char c){return c=='+'||c=='-'?1:c=='*'||c=='/'?2:c=='^'?4:-1;}
 Poly atom(){ws();if(i>=s.size())fail("missing atom");char c=s[i];
  if(c=='+'||c=='-'){++i;auto v=expr(3);return c=='-'?mul(scalar(C(-1)),v):v;}
  if(c=='('){++i;auto v=expr(0);ws();if(s[i++]!=')')fail("missing )");return v;}
  if(isdigit(c)){size_t b=i;while(i<s.size()&&isdigit(s[i]))++i;return scalar(C(cln::cl_float(cln::cl_I(s.substr(b,i-b).c_str()),cln::default_float_format)));}
  if(isalpha(c)){size_t b=i;while(i<s.size()&&(isalnum(s[i])||s[i]=='_'))++i;auto name=s.substr(b,i-b);ws();
   if(i<s.size()&&(s[i]=='['||s[i]=='(')){char open=s[i],close=open=='['?']':')';size_t end=i;int depth=0;do{if(s[end]==open)++depth;else if(s[end]==close)--depth;++end;}while(end<s.size()&&depth);if(depth)fail("unclosed function");auto label=s.substr(b,end-b);auto found=functions.find(label);if(found!=functions.end()){i=end;return scalar(found->second);}++i;std::vector<C> args;while(true){args.push_back(value(expr(0)));ws();if(s[i]==close){++i;break;}if(s[i++]!=',')fail("missing comma");}
    size_t arity=(name=="Num"||name=="Den"||name=="Log"||name=="Zeta")?1:(name=="PolyLog"||name=="Li22"||name=="Li22new")?2:0;
    if(!arity||args.size()!=arity)fail("unsupported function or arity");
    if(name=="PolyLog"||name=="Zeta"){C order=args[0];if(!cln::zerop(cln::imagpart(order))||order!=C(cln::truncate1(cln::realpart(order))))fail("nonintegral order");}
    std::ostringstream canonical;canonical<<name;for(const C& arg:args)canonical<<'|'<<GiNaC::numeric(arg);auto nk=canonical.str();auto cached=numeric_functions.find(nk);if(cached!=numeric_functions.end()){functions[label]=cached->second;return scalar(cached->second);}
    using namespace GiNaC;ex answer;auto num=[&](int k){return numeric(args.at(k));};
    if(name=="Num"){if(args.size()!=1)fail("Num arity");answer=num(0);}
    else if(name=="Den"){if(args.size()!=1)fail("Den arity");answer=1/num(0);}
    else if(name=="Log")answer=log(num(0));
    else if(name=="Zeta")answer=zeta(numeric(cln::truncate1(cln::realpart(args.at(0)))));
    else if(name=="PolyLog")answer=Li(numeric(cln::truncate1(cln::realpart(args.at(0)))),num(1));
    else if((name=="Li22"||name=="Li22new")){lst orders,values;orders.append(2);orders.append(2);values.append(num(0));values.append(num(1));answer=Li(orders,values);}
    else fail("unknown function "+name);
    answer=answer.evalf();if(!is_a<numeric>(answer))fail("unevaluated "+name);C result=ex_to<numeric>(answer).to_cl_N();numeric_functions[nk]=result;functions[label]=result;return scalar(result);
   }
   if(!symbols.count(name))fail("unknown symbol "+name);return symbols.at(name);
  }fail("bad atom");
 }
 Poly expr(int low){auto v=atom();while(true){ws();if(i>=s.size())break;char op=s[i];int p=prec(op);if(p<low)break;++i;auto b=expr(p+(op=='^'?0:1));if(op=='+')v=add(v,b);else if(op=='-')v=add(v,b,-1);else if(op=='*')v=mul(v,b);else if(op=='/')v=divp(v,b);else v=power(v,b);}return v;}
};

int main(int argc,char**argv){
 if(argc!=4)throw std::runtime_error("EXPRESSION SYMBOLS DIGITS");
 int digits=std::stoi(argv[3]);GiNaC::Digits=digits;cln::default_float_format=cln::float_format(digits);
 std::map<std::string,Poly> syms{{"CA",{{{1,0},one()}}},{"CF",{{{0,1},one()}}},{"N",{{{1,0},one()}}}};
 std::ifstream input(argv[2]);std::string name,real,imag;
 while(input>>name>>real>>imag){C a=GiNaC::numeric(real.c_str()).to_cl_N(),b=GiNaC::numeric(imag.c_str()).to_cl_N();syms[name]=scalar(a+C(cln::complex(0,1))*b);}
 std::ifstream in(argv[1]);std::string expression((std::istreambuf_iterator<char>(in)),{});
 Parser p(expression,syms);auto result=p.expr(0);p.ws();if(p.i!=expression.size())throw std::runtime_error("trailing input");
 for(const auto&term:result)std::cout<<term.first.first<<" "<<term.first.second<<" "<<GiNaC::numeric(cln::realpart(term.second))<<" "<<GiNaC::numeric(cln::imagpart(term.second))<<"\n";
}
