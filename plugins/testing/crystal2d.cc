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
 AssignParameter("a", "1");
 AssignParameter("b", "1");
 AssignParameter("gamma", "90");
 AssignParameter("nx", "1");
 AssignParameter("ny", "1");
 AssignParameter("symbol", "H");
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
 std::cout << "      a             : Largo del Vector a                                       \n";
 std::cout << "      b             : Largo del Vector b                                       \n";
 std::cout << "      gamma         : Angulo entre los vectores, en grad.                      \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystalfcc symbol=Ar nx=2 ny=2 nz=2                                     \n";
 std::cout << " input crystalfcc Ar 2 2 2                                                   \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo SC en la simulacion.  \n";
}

void Crystal2DGenerator::Generate(Configuration & conf) const
{
// Vector tmp=sc.GetCell()[2];
 long int cc = 0;
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 std::cerr << "gamma = " << gamma << " a = " << a << " b = " << b << '\n';
 std::cerr << "ax = " << ax << " ay = " << ay << " tan(gamma) = " << tan(gamma) << '\n';
 std::cerr << "sin(gamma) = " << sin(gamma) << " cos(gamma) = " << cos(gamma) << '\n';
 Vector baseA(ax,0,0);
 Vector baseB(ay*cos(gamma),ay*sin(gamma),0);
 std::cerr << "Vector A = " << baseA << '\n';
 std::cerr << "Vector B = " << baseB << '\n';
 for (long j=0;j<ny;++j)
  for (long i=0;i<nx;++i)
  {
//   p = j*baseA + i*baseB;
//   sc.Create(new Atom(spc));
//   sc.SetFracPosition(cc++, p);
   atoms[cc++].Position() = cell.ScaleByCell(j*baseA + i*baseB);
  }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Crystal2DGenerator(args); }
void destroy(Plugin * m) { delete m; }


