//
//
//

#include <lpmd/physunits.h>

#include "pressure.h"

using namespace lpmd;

Pressure::Pressure(std::string args): Module("pressure") 
{ 
 ProcessArguments(args);
}

Pressure::~Pressure() { }

void Pressure::Evaluate(SimulationCell & sc, Potential & pot)
{
 double v = sc.Volume();
 double K=0.0e0, virial = 0.0e0;
 for (long i=0;i<sc.Size();++i)
 {
  Atom a = sc.GetAtom(i);
  double m = a.Mass();
  Vector vel = a.Velocity();
  Vector ff = a.Acceleration()*m;
  virial += (Dot(a.Position(), ff));
  K += m*vel.Mod2();
 }
 K = 0.5*K;
 press = ((1.0/3.0)*virial*KIN2EV + (2.0/3.0)*K)*PRESSFACTOR / v;
}

const double & Pressure::Value() const { return press; }

void Pressure::SetParameter(std::string name) { }

void Pressure::Show() const
{
 Module::Show();
}

std::string Pressure::Keywords() const { return ""; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Pressure(args); }
void destroy(Module * m) { delete m; }


