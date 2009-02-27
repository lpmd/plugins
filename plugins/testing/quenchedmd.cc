//
//
//

#include "quenchedmd.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

QuenchedMDModifier::QuenchedMDModifier(std::string args): Module("quenchedmd")
{
 ProcessArguments(args);
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
}

QuenchedMDModifier::~QuenchedMDModifier() { }

void QuenchedMDModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = quenchedmd                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use quenchedmd                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply quenchedmd start=0 each=1 end=1000                                    \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string QuenchedMDModifier::Keywords() const
{
 return "start end each";
}

void QuenchedMDModifier::Apply(SimulationCell & sc)
{
 //
 for (unsigned long int i=0;i<sc.size();++i)
 {
  if (Dot(sc[i].Acceleration(), sc[i].Velocity()) < 0.0) sc.SetVelocity(i, Vector(0.0, 0.0, 0.0));
 }
}

void QuenchedMDModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 Apply(sc);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new QuenchedMDModifier(args); }
void destroy(Module * m) { delete m; }

