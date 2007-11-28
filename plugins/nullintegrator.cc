//
//
//

#include "nullintegrator.h"

using namespace lpmd;

NullIntegrator::NullIntegrator(std::string args): Module("nullintegrator")
{
 start_step = 1;
 ProcessArguments(args);
}

NullIntegrator::~NullIntegrator() { }

void NullIntegrator::Advance(SimulationCell & sc, Potential & p) { p.UpdateForces(sc); }

void NullIntegrator::SetParameter(std::string name) 
{ 
 if (name == "start")
 {
  AssignParameter("start", GetNextWord());
  start_step = GetInteger("start");
 }
}

void NullIntegrator::Show() const
{
 Module::Show();
 std::cout << "start = " << start_step << '\n';
}

std::string NullIntegrator::Keywords() const { return "start"; }


// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullIntegrator(args); }
void destroy(Module * m) { delete m; }


