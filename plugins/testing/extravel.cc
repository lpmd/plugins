//
//
//

#include "extravel.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

ExtraVelModifier::ExtraVelModifier(std::string args): Module("extravel") 
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 ProcessArguments(args); 
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

ExtraVelModifier::~ExtraVelModifier() { }

void ExtraVelModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para ...                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use extravel                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply extravel start=0 each=10 end=100                                      \n\n";
 std::cout << "      De esta forma aplicamos extravel entre 0 y 100 cada 10 steps.            \n";
}

void ExtraVelModifier::Apply(SimulationCell & sc)
{
 DebugStream() << "-> Applying extravel modifier!" << '\n';
 for (unsigned long int i=0;i<sc.size();++i)
 {
  const Atom & at = sc[i];
  if (at.IsTypeSet() && at.Type().GetBool("extravel"))
  {
   AtomType & atype = at.Type();
   const Vector & vapply = Vector(atype.GetDouble("vx"), atype.GetDouble("vy"), atype.GetDouble("vz"));
   sc.SetVelocity(i, at.Velocity()+vapply);
  }
 }
}

void ExtraVelModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 Apply(sc);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ExtraVelModifier(args); }
void destroy(Module * m) { delete m; }

