//
//
//

#include "varstep.h"

#include <lpmd/simulation.h>
#include <lpmd/session.h>

#include <iostream>

using namespace lpmd;

VariableStep::VariableStep(std::string args): Plugin("varstep", "1.0")
{
 DefineKeyword("dt", "1.0");
 DefineKeyword("scale", "0.15");
 DefineKeyword("start", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 dt = double((*this)["dt"]);
 scalefactor = double((*this)["scale"]);
 start = int((*this)["start"]);
 vhalf = NULL;
}

VariableStep::~VariableStep() 
{ 
 delete vhalf; 
 DebugStream() << "-> Total simulated time: " << elapsed_time << " fs\n";
}

void VariableStep::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para integrar utilizando un metodo de paso variable.\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      dt            : Especifica el tiempo en femto-segundos para el           \n" ;
 std::cout << "                      integrador.                                              \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use varstep                                                                   \n";
 std::cout << "     dt 10.0                                                                   \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " integrator varstep start=1000                                                 \n\n";
 std::cout << "      El integrador puede ser llamado desde el principio (sin usar start) o en \n";
 std::cout << " otro instante de tiempo, para poder modificar el integrador durante la        \n";
 std::cout << " simulacion.                                                                   \n";
}

void VariableStep::Initialize(Simulation & sim, Potential & p) { assert(&p != 0); UseOldConfig(sim); }

void VariableStep::Advance(Simulation & sim, Potential & p)
{
 p.UpdateForces(sim);
 BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicCell & cell = sim.Cell();
 lpmd::BasicParticleSet & oldatoms = OldConfig().Atoms(); 
 const Vector aczero(0.0, 0.0, 0.0);

 long int i = 0;
 Atom at("H");
 // Setea a cero las aceleraciones de los atomos con fixedvel
 if (atoms.HaveAny(Tag("fixedvel"))) 
 {
  GlobalSession.DebugStream() << "-> Considering the fixedvel flag on some atoms\n";
  for (i=0;i<atoms.Size();++i)
  {
   at = atoms[i];
   if (atoms.Have(at, Tag("fixedvel"))) atoms[i].Acceleration() = aczero;
  }
 }

 // Compute vhalf
 if (vhalf == NULL)
 {
  vhalf = new Vector[atoms.Size()];
  elapsed_time = 0.0;
 }
 for (long int i=0;i<atoms.Size();++i)
     vhalf[i] = atoms[i].Velocity() + 0.5*curr_dt*atoms[i].Acceleration();

 // Find the maximal acceleration change
 double max_acc = -1.0, med_acc = 0.0;
 for (long int i=0;i<atoms.Size();++i)
 {
  double dacc = (atoms[i].Acceleration()-oldatoms[i].Acceleration()).Module();
  if (dacc > max_acc) max_acc = dacc;
  med_acc += atoms[i].Acceleration().Module();
 }
 med_acc /= double(atoms.Size());
 
 // do something with the timestep
 curr_dt = scalefactor*(med_acc/max_acc)*dt;
 if (atoms.HaveAny(Tag("fixedpos")))
 {
  GlobalSession.DebugStream() << "-> Considering the fixedpos flag on some atoms\n";
  for (i=0;i<atoms.Size();++i)
  {
   at = atoms[i];
   if (atoms.Have(at, Tag("fixedpos")) && (atoms.GetTag(at, Tag("fixedpos")) == "true")) continue;
   else 
   {
    Vector newpos = atoms[i].Position() + vhalf[i]*curr_dt;
    oldatoms[i].Acceleration() = atoms[i].Acceleration();
    atoms[i].Position() = cell.FittedInside(newpos);
   }
  }
 }
 else
 {
  for (i=0;i<atoms.Size();++i)
  {
   Vector newpos = atoms[i].Position() + vhalf[i]*curr_dt;
   oldatoms[i].Acceleration() = atoms[i].Acceleration();
   atoms[i].Position() = cell.FittedInside(newpos);
  }
 }

 p.UpdateForces(sim);

 // Setea a cero las aceleraciones de los atomos con fixedvel
 if (atoms.HaveAny(Tag("fixedvel"))) 
 {
  GlobalSession.DebugStream() << "-> Considering the fixedvel flag on some atoms\n";
  for (i=0;i<atoms.Size();++i)
  {
   at = atoms[i];
   if (atoms.Have(at, Tag("fixedvel"))) atoms[i].Acceleration() = aczero;
  }
 }

 if (atoms.HaveAny(Tag("fixedpos")))
 {
  GlobalSession.DebugStream() << "-> Considering the fixedpos flag on some atoms\n";
  for (i=0;i<atoms.Size();++i)
  {
   at = atoms[i];
   if (atoms.Have(at, Tag("fixedpos")) && (atoms.GetTag(at, Tag("fixedpos")) == "true")) continue;
   else
     atoms[i].Velocity() = vhalf[i] + 0.5*curr_dt*atoms[i].Acceleration();
  }
 }
 else
 {
  for (i=0;i<atoms.Size();++i)
     atoms[i].Velocity() = vhalf[i] + 0.5*curr_dt*atoms[i].Acceleration();
 }
 elapsed_time += curr_dt;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new VariableStep(args); }
void destroy(Plugin * m) { delete m; }

