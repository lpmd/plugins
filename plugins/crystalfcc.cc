//
//
//

#include "crystalfcc.h"

using namespace lpmd;

FCCGenerator::FCCGenerator(std::string args): Module("crystalfcc")
{
 nx = 1;
 ny = 1;
 nz = 1;
 spc = 1;
 a = 1.0;
 ProcessArguments(args);
}

FCCGenerator::~FCCGenerator() { }

void FCCGenerator::SetParameter(std::string name)
{
 if (name == "symbol") 
 {
  AssignParameter("symbol", GetNextWord());
  spc = ElemNum(GetString("symbol"));
 }
 else if (name == "nx") 
 {
  AssignParameter("nx", GetNextWord());
  nx = GetInteger("nx");
 }
 else if (name == "ny")
 {
  AssignParameter("ny", GetNextWord());
  ny = GetInteger("ny");
 }
 else if (name == "nz")
 {
  AssignParameter("nz", GetNextWord());
  nz = GetInteger("nz");
 }
 else if (name == "a")
 {
  AssignParameter("a", GetNextWord());
  a = GetDouble("a");
 }
}

void FCCGenerator::Show() const
{
 Module::Show();
 std::cout << "   symbol = " << ElemSym[spc] << '\n';
 std::cout << "       nx = " << nx << '\n';
 std::cout << "       ny = " << ny << '\n';
 std::cout << "       nz = " << nz << '\n';
 std::cout << "        a = " << a << '\n';
}

std::string FCCGenerator::Keywords() const
{
 return "symbol nx ny nz a";
}

void FCCGenerator::Generate(SimulationCell & sc) const
{
 Vector p;
 long int cc = 0;
 for (long k=0;k<nz;++k)
 {
  for (long j=0;j<ny;++j)
  {
   for (long i=0;i<nx;++i)
   {
    p = Vector((double(i)+0.5)*a, (double(j)+0.5)*a, (double(k)+0.5)*a);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+0.5)*a, double(j)*a, (double(k)+1.0)*a);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+1.0)*a, double(j)*a, (double(k)+0.5)*a);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
    p = Vector((double(i)+1.0)*a, (double(j)+0.5)*a, (double(k)+1.0)*a);
    sc.AppendAtom(Atom(spc));
    sc.SetFracPosition(cc++, p);
   }
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new FCCGenerator(args); }
void destroy(Module * m) { delete m; }


