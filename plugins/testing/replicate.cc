//
//
//

#include "replicate.h"

#include <lpmd/simulation.h>

#include <iostream>

using namespace lpmd;

ReplicateModifier::ReplicateModifier(std::string args): Module("replicate")
{
 ParamList & param = (*this);
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("nx", "1");
 DefineKeyword("ny", "1");
 DefineKeyword("nz", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 nx = int(param["nx"]);
 ny = int(param["ny"]);
 nz = int(param["nz"]);
}

ReplicateModifier::~ReplicateModifier() { }

void ReplicateModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << "  Aplicando el Modulo                                                          \n";
 std::cout << " prepare replicate nx=3 ny=3 nz=3                                              \n";
}

void ReplicateModifier::Apply(Simulation & con)
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 unsigned long int Ntmp = atoms.Size();
 //Atom * atomos;
 //atomos = new Atom[Ntmp];
 //for (unsigned long int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 //for (unsigned long int i=0;i<Ntmp;i++) { sc.Create(new Atom(atomos[i]));}
 for (unsigned long i=1;i<nx;i++)
 {
  for (unsigned long int j=0;j<Ntmp;j++)
  {
   std::string symb = atoms[j].Symbol();
   lpmd::Vector pos = atoms[j].Position() + cell[0]*i;
   lpmd::Vector vel = atoms[j].Velocity();
   lpmd::Vector acc = atoms[j].Acceleration();
   lpmd::Atom atm = Atom(symb,pos,vel,acc);
   atoms.Append(atm);
  }
 }
 Ntmp = atoms.Size();
 for (unsigned long i=1;i<ny;i++)
 {
  for(unsigned long int j=0;j<Ntmp;j++)
  {
   std::string symb = atoms[j].Symbol();
   lpmd::Vector pos = atoms[j].Position() + cell[0]*i;
   lpmd::Vector vel = atoms[j].Velocity();
   lpmd::Vector acc = atoms[j].Acceleration();
   lpmd::Atom atm = Atom(symb,pos,vel,acc);
   atoms.Append(atm);
  }
 }
 Ntmp = atoms.Size();
 for (unsigned long i=1;i<nz;i++)
 {
  for (unsigned long int j=0;j<Ntmp;j++)
  {
   std::string symb = atoms[j].Symbol();
   lpmd::Vector pos = atoms[j].Position() + cell[0]*i;
   lpmd::Vector vel = atoms[j].Velocity();
   lpmd::Vector acc = atoms[j].Acceleration();
   lpmd::Atom atm = Atom(symb,pos,vel,acc);
   atoms.Append(atm);
  }
 }
 //Resetea los vectores base de la celda.
 cell[0] = cell[0]*nx;
 cell[1] = cell[1]*ny;
 cell[2] = cell[2]*nz;
 //Asigna el index() a cada atomo de la celda.
#warning no hay clearforces ni assignindex
 //sc.AssignIndex();
 //sc.ClearForces();
}


// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ReplicateModifier(args); }
void destroy(Module * m) { delete m; }


