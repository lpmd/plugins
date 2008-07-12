//
//
//

#include "moleculecm.h"

#include <lpmd/md.h>
#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

MoleculeCMModifier::MoleculeCMModifier(std::string args): Module("moleculecm")
{
 AssignParameter("radius", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 radius = GetDouble("radius");
}

MoleculeCMModifier::~MoleculeCMModifier() { }

void MoleculeCMModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = moleculecm                                               \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use moleculecm                                                                \n";
 std::cout << "     radius 1.1                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply moleculecm                                                            \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string MoleculeCMModifier::Keywords() const
{
 return "radius";
}

void MoleculeCMModifier::Apply(SimulationCell & sc)
{
 std::list<Neighbor> nlist;
 std::list<Atom> tmplist;
 int * used = new int[sc.Size()];
 for (int i=0;i<sc.Size();++i) used[i] = 0;
 for (long i=0;i<sc.Size();++i) 
 {
  nlist.clear();
  Vector cmpos; 
  double m = 0.0;
  if (used[i] == 0)
  {
   cmpos += (sc[i].Mass()*sc[i].Position());
   m += sc[i].Mass();
   sc.BuildNeighborList(i, nlist, true, radius);
   Neighbor * closest = NULL;
   for (std::list<Neighbor>::iterator it=nlist.begin();it!=nlist.end();++it)
   {
    Neighbor & nn = *it;
    if (nn.r < radius)
    {
     if (used[nn.j->Index()] == 0)
     {
      if (closest == NULL) closest = &nn;
      else if (nn.r < closest->r) closest = &nn;
     }
    }
   }
   if (closest != NULL)
   {
    cmpos = cmpos + ((closest->j)->Mass()*(sc[i].Position()+closest->rij));
    m += (closest->j)->Mass();
    cmpos = (1.0/m)*cmpos;
    tmplist.push_back(Atom(sc[i].Species(), cmpos));
    used[i] = 1;
    used[closest->j->Index()] = 1;
   }
  }
 }
 sc.Clear();
 for (std::list<Atom>::const_iterator it=tmplist.begin();it!=tmplist.end();++it) sc.AppendAtom(*it);
 delete [] used;
}

void MoleculeCMModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MoleculeCMModifier(args); }
void destroy(Module * m) { delete m; }


