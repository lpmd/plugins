//
//
//

#include "metropolis.h"

#include <lpmd/simulation.h>
#include <lpmd/properties.h>
#include <lpmd/session.h>

#include <cmath>
#include <iostream>

using namespace lpmd;

Metropolis::Metropolis(std::string args): Plugin("metropoli", "2.0")
{
 DefineKeyword("start", "1");
 DefineKeyword("temp","0.0");
 DefineKeyword("percent","5");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int((*this)["start"]);
 Temp = double((*this)["temp"]);
 percent = double((*this)["percent"]);
}

Metropolis::~Metropolis() { }

void Metropolis::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para minimzar utilizando el metodo de metropoli.  \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      temp            Utilizado para asignar una temperatura externa.          \n";
 std::cout << "      rand            Utilizado para el tamaño de los vectores al azar.        \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use metropoli                                                                 \n";
 std::cout << "    temp 300                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << '\n';
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator metropoli start=1 end=500 each=2                                 \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void Metropolis::Initialize(Simulation & sim, Potential & p) { UseOldConfig(sim); }

void Metropolis::Advance(Simulation & sim, long i)
{
 const double kboltzmann = double(GlobalSession["kboltzmann"]);
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::CombinedPotential & pots = sim.Potentials();
 lpmd::BasicCell & cell = sim.Cell();
 double penergy = 0.0e0;
 double min = sim.MinimumPairDistance();
 for(int j=0;j<pots.Size();++j)
 {
  penergy += pots[j].energy(sim);
 }
 Vector oldpos = atoms[i].Position();
 Vector random = RandomVector(min*percent/100);
 Vector newpos = cell.FittedInside(atoms[i].Position() + random);
 atoms[i].Position() = newpos;
 pots.UpdateForces(sim);
 double npe = 0.0e0;
 for(int j=0;j<pots.Size();++j)
 {
  npe += pots[j].energy(sim);
 }
 if ((npe)<(penergy))
 {
  atoms[i].Position() = newpos;
 }
 else if(exp(-(npe-penergy)/kboltzmann*Temp)<1)
 {
  atoms[i].Position() = newpos;
 }
 else
 {
  atoms[i].Position() = oldpos;
 }
 pots.UpdateForces(sim);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Metropolis(args); }
void destroy(Plugin * m) { delete m; }
