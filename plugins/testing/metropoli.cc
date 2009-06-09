//
//
//

#include "metropoli.h"

#include <lpmd/simulation.h>
#include <lpmd/properties.h>

#include <cmath>
#include <iostream>

using namespace lpmd;

Metropoli::Metropoli(std::string args): Plugin("metropoli", "2.0")
{
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int((*this)["start"]);
 Temp = double((*this)["temp"]);
}

Metropoli::~Metropoli() { }

void Metropoli::ShowHelp() const
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

void Metropoli::Initialize(Simulation & sim, Potential & p) { UseOldConfig(sim); }

void Metropoli::Advance(Simulation & sim, long i)
{
 const double kboltzmann = double(Parameter(sim.GetTag(sim,"kboltzmann")));
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 //lpmd::BasicParticleSet & oldatoms = OldConfig().Atoms();
 lpmd::CombinedPotential & pots = sim.Potentials();
 double kenergy = KineticEnergy(atoms);
 double penergy = 0.0e0;
 double min = sim.MinimumPairDistance();
 for(int j=0;j<pots.Size();++j)
 {
  penergy += pots[j].energy(sim);
 }
 Vector oldpos = atoms[i].Position();
 std::cerr << "min = " << min*20/100 << '\n';
 Vector random = RandomVector(min*20/100);
 Vector newpos = atoms[i].Position() + random;
 atoms[i].Position() = newpos;
 double nke = KineticEnergy(atoms);
 double npe = 0.0e0;
 for(int j=0;j<pots.Size();++j)
 {
  npe += pots[j].energy(sim);
 }
 if ((nke+npe)<(kenergy+penergy))
 {
  atoms[i].Position() = newpos;
 }
 else if(exp(-(nke+npe-kenergy-penergy)/kboltzmann*Temp)<1)
 {
  atoms[i].Position() = newpos;
 }
 else
 {
  atoms[i].Position() = oldpos;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Metropoli(args); }
void destroy(Plugin * m) { delete m; }
