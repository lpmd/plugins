//
//
//

#include "displace.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

DisplaceModifier::DisplaceModifier(std::string args): Module("displace")
{
 AssignParameter("x", "1.0");
 AssignParameter("y", "0.0");
 AssignParameter("z", "0.0");
 AssignParameter("autocenter", "false");
 // 
 ProcessArguments(args);
 offset = Vector(GetDouble("x"), GetDouble("y"), GetDouble("z"));
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

DisplaceModifier::~DisplaceModifier() { }

void DisplaceModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = displace                                                 \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      x             : Coord. X del vector de desplazamiento                    \n";
 std::cout << "      y             : Coord. Y del vector de desplazamiento                    \n";
 std::cout << "      z             : Coord. Z del vector de desplazamiento                    \n";
 std::cout << "      autocenter    : Si es true se escoge x, y, z automaticamente para centrar\n";
 std::cout << "                      los atomos en la celda de simulacion (default=false)     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " prepare displace x=1.0 y=0.0 z=0.0                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string DisplaceModifier::Keywords() const
{
 return "x y z autocenter start end each";
}

void DisplaceModifier::Apply(SimulationCell & sc)
{
 if (GetString("autocenter") == "true")
 {
  Vector cm(0.0, 0.0, 0.0);
  Vector boxcenter(0.5, 0.5, 0.5);
  sc.ScaleByCell(boxcenter);
  double tm = 0.0;
  for (long int i=0;i<sc.Size();++i)
  {
   const Atom & at = sc[i];
   cm += at.Mass()*at.Position();
   tm += at.Mass();
  }
  cm = cm*(1.0/tm);
  offset = boxcenter - cm;
  std::cerr << "-> Autocentering atoms: offset = " << offset << '\n';
 }
 for (long int i=0;i<sc.Size();++i)
 {
  Vector pos = sc.GetAtom(i).Position() + offset; 
  sc.SetPosition(i, pos);
 }
}

void DisplaceModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 Apply(sc);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new DisplaceModifier(args); }
void destroy(Module * m) { delete m; }

