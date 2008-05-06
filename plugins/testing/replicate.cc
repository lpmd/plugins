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
 AssignParameter("nx", "1");
 AssignParameter("ny", "1");
 AssignParameter("nz", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 nx = GetInteger("nx");
 ny = GetInteger("ny");
 nz = GetInteger("nz");
}

ReplicateModifier::~ReplicateModifier() { }

void ReplicateModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = replicate                                                \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << "  Aplicando el Modulo                                                          \n";
 std::cout << " prepare replicate nx=3 ny=3 nz=3                                              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string ReplicateModifier::Keywords() const
{
 return "nx ny nz";
}

void ReplicateModifier::Apply(SimulationCell & sc)
{
 int Ntmp = sc.Size();
 Atom *atomos;
 atomos = new Atom[Ntmp];
 for(int i=0;i<Ntmp;i++) atomos[i]=sc.GetAtom(i);
 sc.Initialize(nx*Ntmp);
 //int counter = 0;
 for(int i=0;i<Ntmp;i++){sc.AppendAtom(atomos[i]);}//counter++;}

 for(unsigned long i=1;i<nx;i++)
 {
  for(int j=0;j<Ntmp;j++)
  {
   Atom tmp=atomos[j];
   tmp.SetPos(atomos[j].Position()+sc.GetVector(0)*i);
   sc.AppendAtom(tmp);
  }
 }
 delete[] atomos;
 Ntmp = sc.Size();
 atomos = new Atom[Ntmp];
 for(int i=0;i<Ntmp;i++) atomos[i]=sc.GetAtom(i);
 sc.Initialize(ny*Ntmp);
 for(int i=0;i<Ntmp;i++){sc.AppendAtom(atomos[i]);}
 for(unsigned long i=1;i<ny;i++)
 {
  for(int j=0;j<Ntmp;j++)
  {
   Atom tmp=atomos[j];
   tmp.SetPos(atomos[j].Position()+sc.GetVector(1)*i);
   sc.AppendAtom(tmp);
  }
 }
 delete[] atomos;
 Ntmp = sc.Size();
 atomos = new Atom[Ntmp];
 for(int i=0;i<Ntmp;i++) atomos[i]=sc.GetAtom(i);
 sc.Initialize(nz*Ntmp);
 for(int i=0;i<Ntmp;i++) { sc.AppendAtom(atomos[i]);}
 for(unsigned long i=1;i<nz;i++)
 {
  for(int j=0;j<Ntmp;j++)
  {
   Atom tmp=atomos[j];
   tmp.SetPos(atomos[j].Position()+sc.GetVector(2)*i);
   sc.AppendAtom(tmp);
  }
 }
 delete[] atomos;
 //Resetea los vectores base de la celda.
 Vector a=sc.GetVector(0);
 sc.SetVector(0, a*nx);
 Vector b=sc.GetVector(1);
 sc.SetVector(1, b*ny);
 Vector c=sc.GetVector(2);
 sc.SetVector(2, c*nz);
 //Asigna el index() a cada atomo de la celda.
 sc.AssignIndex();
 sc.ClearForces();
}

void ReplicateModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ReplicateModifier(args); }
void destroy(Module * m) { delete m; }


