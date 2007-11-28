//
//
//

#include <cmath>

#include "skewstart.h"

using namespace lpmd;

SkewStartGenerator::SkewStartGenerator(std::string args): Module("skewstart")
{
 ProcessArguments(args);
}

SkewStartGenerator::~SkewStartGenerator() { }

void SkewStartGenerator::SetParameter(std::string name)
{
 if (name == "symbol") 
 {
  AssignParameter("symbol", GetNextWord());
  spc = ElemNum(GetString("symbol"));
 }
 else if (name == "atoms")
 {
  AssignParameter("atoms", GetNextWord());
  n = GetInteger("atoms");
 }
}

void SkewStartGenerator::Show() const
{
 Module::Show();
 std::cout << "   symbol = " << ElemSym[spc] << '\n';
 std::cout << "   atoms = " << n << '\n';
}

std::string SkewStartGenerator::Keywords() const
{
 return "symbol atoms";
}

//
// Skew Start method, like Moldy
// Only works with atoms (we have no molecules yet!)
//
void SkewStartGenerator::Generate(SimulationCell & sc) const
{
 int h, k, l;
 double dx, dy, dz;
 h = int(pow(double(n), 2.0/3.0));
 k = int(pow(double(n), 1.0/3.0));
 l = 1;
 dx = h / double(n);
 dy = k / double(n);
 dz = l / double(n);
 for (long i=0;i<n;++i)
 {
  Atom a(spc, Vector(0.0));
  sc.AppendAtom(a);
  sc.SetFracPosition(i, Vector(dx*double(i), dy*double(i), dz*double(i)));
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new SkewStartGenerator(args); }
void destroy(Module * m) { delete m; }




