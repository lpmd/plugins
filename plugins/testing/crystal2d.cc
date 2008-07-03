//
//
//

#include "crystal2d.h"

#include <lpmd/atom.h>
#include <lpmd/simulationcell.h>

using namespace lpmd;

Crystal2DGenerator::Crystal2DGenerator(std::string args): Module("crystal2d")
{
 AssignParameter("a", "1");
 AssignParameter("b", "1");
 AssignParameter("gamma", "90");
 AssignParameter("nx", "1");
 AssignParameter("ny", "1");
 AssignParameter("symbol", "H");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 a = GetDouble("a");
 b = GetDouble("b");
 gamma = M_PI*GetDouble("gamma")/180.0e0;
 spc = ElemNum(GetString("symbol"));
 nx = GetInteger("nx");
 ny = GetInteger("ny");
 nz = GetInteger("nz");
}

Crystal2DGenerator::~Crystal2DGenerator() { }

void Crystal2DGenerator::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = crystal2d                                                \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para crear celdas bidimensionales.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      a             : Largo del Vector a                                       \n";
 std::cout << "      b             : Largo del Vector b                                       \n";
 std::cout << "      gamma         : Angulo entre los vectores, en grad.                      \n";
 std::cout << "      symbol        : Especifica la especie atomica, utilizando su simbolo.    \n";
 std::cout << "      nx            : Repeticiones en la direccion X.                          \n";
 std::cout << "      ny            : Repeticiones en la direccion Y.                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Utilizando el Modulo :                                                        \n";
 std::cout << " input crystalfcc symbol=Ar nx=2 ny=2 nz=2                                     \n";
 std::cout << " input crystalfcc Ar 2 2 2                                                   \n\n";
 std::cout << "      De esta forma creamos una celda de entrada de tipo SC en la simulacion.  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Crystal2DGenerator::Keywords() const
{
 return "a b gamma symbol nx ny";
}

void Crystal2DGenerator::Generate(SimulationCell & sc) const
{
 Vector tmp=sc.GetVector(2);
 Vector p;
 long int cc = 0;
 double ax = 1.0/double(nx);
 double ay = 1.0/double(ny);
 std::cerr << "gamma = " << gamma << " a = " << a << " b = " << b << '\n';
 std::cerr << "ax = " << ax << " ay = " << ay << " tan(gamma) = " << tan(gamma) << '\n';
 std::cerr << "sin(gamma) = " << sin(gamma) << " cos(gamma) = " << cos(gamma) << '\n';
 Vector baseA(a,0,0);
 Vector baseB(b*cos(gamma),b*sin(gamma),0);
 std::cerr << "Vector A = " << baseA << '\n';
 std::cerr << "Vector B = " << baseB << '\n';
 for (long j=0;j<ny;++j)
 {
  for (long i=0;i<nx;++i)
  {
//   if((j%2)!=0) p = Vector((double(i)*ax+double(j)*ay/tan(gamma)), (double(j))*ay, 0.0e0);
//   else p = Vector((double(i)*ax),(double(j)*ay), 0.0e0);
   p = Vector(double(i)*ax,double(j)*ay,0);
   sc.AppendAtom(Atom(spc));
   sc.SetFracPosition(cc++, p);
   p = Vector(double(i)*ax+(double(j)+0.5)*ay/tan(gamma),(double(j)+0.5)*ay,0);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Crystal2DGenerator(args); }
void destroy(Module * m) { delete m; }


