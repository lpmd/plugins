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
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("radius", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 radius = GetDouble("radius");
}

MoleculeCMModifier::~MoleculeCMModifier() { }

void MoleculeCMModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use moleculecm                                                                \n";
 std::cout << "     radius 1.1                                                                \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply moleculecm                                                            \n\n";
}

void MoleculeCMModifier::Apply(SimulationCell & sc)
{
 std::vector<Neighbor> nlist;
 std::list<Atom> tmplist;
 int * used = new int[sc.size()];
 for (unsigned long int i=0;i<sc.size();++i) used[i] = 0;
 for (unsigned long int i=0;i<sc.size();++i) 
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
   for (unsigned long int k=0;k<nlist.size();++k)
   {
    Neighbor & nn = nlist[k];
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
 sc.clear(); // FIXME debe borrar los atomos en memoria tambien 
 for (std::list<Atom>::const_iterator it=tmplist.begin();it!=tmplist.end();++it) sc.Create(new Atom(*it));
 delete [] used;
}

void MoleculeCMModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MoleculeCMModifier(args); }
void destroy(Module * m) { delete m; }


