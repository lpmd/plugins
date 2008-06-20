//
//
//

#include "rotate.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

RotateModifier::RotateModifier(std::string args): Module("rotate")
{
 AssignParameter("x", "1.0");
 AssignParameter("y", "0.0");
 AssignParameter("z", "0.0");
 AssignParameter("angle", "0.0");
 // 
 ProcessArguments(args);
 axis = Vector(GetDouble("x"), GetDouble("y"), GetDouble("z"));
 axis.Norm();
 angle = M_PI*GetDouble("angle")/180.0;
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 //
 double c=cos(angle), s=sin(angle), t=1.0-c;
 double x=axis.Get(0), y=axis.Get(1), z=axis.Get(2);
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = rotate                                                   \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para rotar el sistema en torno a un eje en un     \n";
 std::cout << " angulo determinado. Puede aplicarse al inicio en la instruccion \"prepare\"   \n";
 std::cout << " como durante la simulacion en la instruccion \"apply\".                       \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      x             : Coord. X del eje de rotacion                             \n";
 std::cout << "      y             : Coord. Y del eje de rotacion                             \n";
 std::cout << "      z             : Coord. Z del eje de rotacion                             \n";
 std::cout << "      angle         : Angulo de rotacion (en grados)                           \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare rotate x=1.0 y=0.0 z=0.0 angle=45.0                                   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string RotateModifier::Keywords() const
{
 return "x y z angle start end each";
}

void RotateModifier::Apply(SimulationCell & sc)
{
 Vector center = (sc.GetVector(0)+sc.GetVector(1)+sc.GetVector(2))*0.5;
 for (long int i=0;i<sc.Size();++i)
 {
  // translate so that the center of the cell is (0, 0, 0)
  Vector pos = sc.GetAtom(i).Position() - center; 
  // rotate
  double rv[3];
  for (int j=0;j<3;++j)
  {
   rv[j] = 0.0;
   for (int i=0;i<3;++i) rv[j] += rotmat[j][i]*pos.Get(i);
   pos.Set(j, rv[j]);
  }
  // untranslate
  sc.SetPosition(i, pos+center);
 }
}

void RotateModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 Apply(sc);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new RotateModifier(args); }
void destroy(Module * m) { delete m; }

