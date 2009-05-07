//
//
//

#include "replicate.h"

#include <lpmd/md.h>
#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

ReplicateModifier::ReplicateModifier(std::string args): Module("replicate")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("nx", "1");
 DefineKeyword("ny", "1");
 DefineKeyword("nz", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 nx = GetInteger("nx");
 ny = GetInteger("ny");
 nz = GetInteger("nz");
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

void ReplicateModifier::Apply(SimulationCell & sc)
{
 unsigned long int Ntmp = sc.size();
 Atom * atomos;
 atomos = new Atom[Ntmp];
 for (unsigned long int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 for (unsigned long int i=0;i<Ntmp;i++) { sc.Create(new Atom(atomos[i]));}
 for (unsigned long i=1;i<nx;i++)
 {
  for (unsigned long int j=0;j<Ntmp;j++)
  {
   Atom * tmp = new Atom(atomos[j]);
   tmp->SetPos(atomos[j].Position()+sc.GetCell()[0]*i);
   sc.Create(tmp);
  }
 }
 delete[] atomos;
 Ntmp = sc.size();
 atomos = new Atom[Ntmp];
 for (unsigned long int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 for (unsigned long int i=0;i<Ntmp;i++) { sc.Create(new Atom(atomos[i]));}
 for (unsigned long i=1;i<ny;i++)
 {
  for(unsigned long int j=0;j<Ntmp;j++)
  {
   Atom * tmp = new Atom(atomos[j]);
   tmp->SetPos(atomos[j].Position()+sc.GetCell()[1]*i);
   sc.Create(tmp);
  }
 }
 delete[] atomos;
 Ntmp = sc.size();
 atomos = new Atom[Ntmp];
 for (unsigned long int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 for (unsigned long int i=0;i<Ntmp;i++) { sc.Create(new Atom(atomos[i]));}
 for (unsigned long i=1;i<nz;i++)
 {
  for (unsigned long int j=0;j<Ntmp;j++)
  {
   Atom * tmp = new Atom(atomos[j]);
   tmp->SetPos(atomos[j].Position()+sc.GetCell()[2]*i);
   sc.Create(tmp);
  }
 }
 delete[] atomos;
 //Resetea los vectores base de la celda.
 Vector a=sc.GetCell()[0];
 sc.GetCell()[0] = a*nx;
 Vector b=sc.GetCell()[1];
 sc.GetCell()[1] = b*ny;
 Vector c=sc.GetCell()[2];
 sc.GetCell()[2] = c*nz;
 //Asigna el index() a cada atomo de la celda.
 sc.AssignIndex();
 sc.ClearForces();
}

void ReplicateModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ReplicateModifier(args); }
void destroy(Module * m) { delete m; }


