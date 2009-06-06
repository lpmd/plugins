//
//
//

#include "skewstart.h"

#include <lpmd/atom.h>
#include <lpmd/configuration.h>
#include <cmath>

using namespace lpmd;

SkewStartGenerator::SkewStartGenerator(std::string args): Plugin("skewstart", "2.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("atoms");
 DefineKeyword("symbol", "H");
 ProcessArguments(args); 
 spc = ElemNum(params["symbol"]);
 n = int(params["atoms"]);
}

SkewStartGenerator::~SkewStartGenerator() { }

void SkewStartGenerator::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para generar una celda inicial con el metodo      \n";
 std::cout << " skewstart disenado por Keith Refson para moldy.                               \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      atoms         : Especifica el numero de atomos para la celda a generar   \n";
 std::cout << "      symbol        : Especifica el simbolo atomico de la especie a generar.   \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=skewstart atoms=108 symbol=Ar                                    \n";
 std::cout << "      De esta forma podemos generar una celda inicial con el metodo skewstart  \n";
 std::cout << " de 108 atomos de Ar.                                                          \n";
}

//
// Skew Start method, like Moldy
// Only works with atoms (we have no molecules yet!)
//
void SkewStartGenerator::Generate(Configuration & config) const
{
 BasicParticleSet & atoms = config.Atoms();
 BasicCell & cell = config.Cell();  
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
//  sc.Create(new Atom(spc));
//  sc.SetFracPosition(i, Vector(dx*double(i), dy*double(i), dz*double(i)));
  atoms[i].Position() = cell.ScaleByCell(Vector(dx*double(i), dy*double(i), dz*double(i)));
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SkewStartGenerator(args); }
void destroy(Plugin * m) { delete m; }
