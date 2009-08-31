//
//
//

#include "rotate.h"

#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

RotateModifier::RotateModifier(std::string args): Plugin("rotate", "2.0")
{
 ParamList & params = (*this);
 DefineKeyword("x", "1.0");
 DefineKeyword("y", "0.0");
 DefineKeyword("z", "0.0");
 DefineKeyword("angle", "0.0");
 // 
 ProcessArguments(args);
 axis = Vector(double(params["x"]), double(params["y"]), double(params["z"]));
 axis.Normalize();
 angle = M_PI*double(params["angle"])/180.0;
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 //
 double c=cos(angle), s=sin(angle), t=1.0-c;
 double x=axis[0], y=axis[1], z=axis[2];
 rotmat[0][0] = t*x*x+c;
 rotmat[0][1] = t*x*y+s*z;
 rotmat[0][2] = t*x*z-s*y;
 rotmat[1][0] = t*x*y-s*z;
 rotmat[1][1] = t*y*y+c;
 rotmat[1][2] = t*y*z+s*x;
 rotmat[2][0] = t*x*z+s*y;
 rotmat[2][1] = t*y*z-s*x;
 rotmat[2][2] = t*z*z+c;
}

RotateModifier::~RotateModifier() { }

void RotateModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para rotar el sistema en torno a un eje en un     \n";
 std::cout << " angulo determinado. Puede aplicarse al inicio en la instruccion \"prepare\"   \n";
 std::cout << " como durante la simulacion en la instruccion \"apply\".                       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      x             : Coord. X del eje de rotacion                             \n";
 std::cout << "      y             : Coord. Y del eje de rotacion                             \n";
 std::cout << "      z             : Coord. Z del eje de rotacion                             \n";
 std::cout << "      angle         : Angulo de rotacion (en grados)                           \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare rotate x=1.0 y=0.0 z=0.0 angle=45.0                                   \n";
}

void RotateModifier::Apply(Configuration & conf)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 Vector center = (cell[0]+cell[1]+cell[2])*0.5;
 for (long i=0;i<atoms.Size();++i)
 {
  // translate so that the center of the cell is (0, 0, 0)
  Vector pos = atoms[i].Position() - center; 
  // rotate
  double rv[3];
  for (int j=0;j<3;++j)
  {
   rv[j] = 0.0;
   for (int k=0;k<3;++k) rv[j] += rotmat[j][k]*pos[k];
   pos[j] = rv[j];
  }
  atoms[i].Position() = cell.FittedInside(pos+center);
 }
}

void RotateModifier::Apply(Simulation & md)
{
 Configuration & sc = md;
 Apply(sc);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new RotateModifier(args); }
void destroy(Plugin * m) { delete m; }

