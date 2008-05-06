//
//
//

#include "skewstart.h"

#include <lpmd/atom.h>
#include <lpmd/simulationcell.h>

#include <cmath>

using namespace lpmd;

SkewStartGenerator::SkewStartGenerator(std::string args): Module("skewstart") 
{ 
 ProcessArguments(args); 
 spc = ElemNum(GetString("symbol"));
 n = GetInteger("atoms");
}

SkewStartGenerator::~SkewStartGenerator() { }

void SkewStartGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = skewstart                                                \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para generar una celda inicial con el metodo      \n";
 std::cout << " skewstart disenado por Keith Refson para moldy.                               \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      atoms         : Especifica el numero de atomos para la celda a generar   \n";
 std::cout << "      symbol        : Especifica el simbolo atomico de la especie a generar.   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=skewstart atoms=108 symbol=Ar                                    \n";
 std::cout << "      De esta forma podemos generar una celda inicial con el metodo skewstart  \n";
 std::cout << " de 108 atomos de Ar.                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
  Atom a(spc, Vector());
  sc.AppendAtom(a);
  sc.SetFracPosition(i, Vector(dx*double(i), dy*double(i), dz*double(i)));
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new SkewStartGenerator(args); }
void destroy(Module * m) { delete m; }
