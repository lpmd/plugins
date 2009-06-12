//
//
//

#include "crystal2d.h"

#include <lpmd/atom.h>
#include <lpmd/configuration.h>

using namespace lpmd;

Crystal2DGenerator::Crystal2DGenerator(std::string args): Plugin("crystal2d", "2.0")
{
 ParamList & params = (*this);
 DefineKeyword("a", "1");
 DefineKeyword("b", "1");
 DefineKeyword("gamma", "90");
 DefineKeyword("nx", "1");
 DefineKeyword("ny", "1");
 DefineKeyword("symbol", "H");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 a = params["a"];
 b = params["b"];
 gamma = M_PI*double(params["gamma"])/180.0e0;
 spc = ElemNum(params["symbol"]);
 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);
}

Crystal2DGenerator::~Crystal2DGenerator() { }

void Crystal2DGenerator::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para crear celdas bidimensionales.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      a             : Largo del vector a                                       \n";
 std::cout << "      b             : Largo del vector b                                       \n";
 std::cout << "      gamma         : Angulo entre los vectores, en grados                     \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystal2d a=1.0 b=1.0 gamma=90.0 symbol=Ar nx=2 ny=2                    \n";
 std::cout << " De esta forma creamos una red cuadrada en la simulacion.                      \n";
}

void Crystal2DGenerator::Generate(Configuration & conf) const
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 assert(atoms.Size() == 0);
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 Vector baseA(ax,0,0);
 Vector baseB(ay*cos(gamma),ay*sin(gamma),0);
 DebugStream() << "-> Using vector a = ";
 FormattedWrite(DebugStream(), baseA);
 DebugStream() << '\n';
 DebugStream() << "-> Using vector b = ";
 FormattedWrite(DebugStream(), baseB);
 DebugStream() << '\n';
 long int cc = 0;
 for (long j=0;j<ny;++j)
  for (long i=0;i<nx;++i)
  {
   atoms.Append(Atom(spc));
   atoms[cc++].Position() = cell.ScaleByCell(j*baseA + i*baseB);
  }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Crystal2DGenerator(args); }
void destroy(Plugin * m) { delete m; }

