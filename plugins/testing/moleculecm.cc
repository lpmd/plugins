//
//
//

#include "moleculecm.h"

#include <lpmd/simulation.h>

#include <iostream>

using namespace lpmd;

MoleculeCMModifier::MoleculeCMModifier(std::string args): Plugin("moleculecm", "2.0")
{
 //
 DefineKeyword("radius", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 radius = double((*this)["radius"]);
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

void MoleculeCMModifier::Apply(Simulation & con)
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();

 std::list<lpmd::Atom> tmplist;
 int * used = new int[atoms.Size()];
 for (long int i=0;i<atoms.Size();++i) used[i] = 0;
 for (long int i=0;i<atoms.Size();++i) 
 {
  //nlist.clear();
  Vector cmpos; 
  double m = 0.0;
  if (used[i] == 0)
  {
   cmpos = cmpos + (atoms[i].Mass()*atoms[i].Position());
   m += atoms[i].Mass();
   lpmd::NeighborList & nlist = con.Neighbors(i,true,radius);
   lpmd::AtomPair * closest = NULL;
   for (long int k=0;k<nlist.Size();++k)
   {
    lpmd::AtomPair & nn = nlist[k];
    if (nn.r < radius)
    {
#warning Index no va más??
//     if (used[nn.j->Index()] == 0)
     {
      if (closest == NULL) closest = &nn;
      else if (nn.r < closest->r) closest = &nn;
     }
    }
   }
   if (closest != NULL)
   {
    cmpos = cmpos + ((closest->j)->Mass()*(atoms[i].Position()+closest->rij));
    m += (closest->j)->Mass();
    cmpos = (1.0/m)*cmpos;
    tmplist.push_back(Atom(atoms[i].Symbol(), cmpos));
    used[i] = 1;
#warning otro Index()
//    used[closest->j->Index()] = 1;
   }
  }
 }
#warning antiguo fixme de borrar atomos en memoria también.
 atoms.Clear(); // FIXME debe borrar los atomos en memoria tambien 
#warning que era esto!! .. no debería actuar ahora directo sobre atoms.?¿?
// for (std::list<Atom>::const_iterator it=tmplist.begin();it!=tmplist.end();++it) sc.Create(new Atom(*it));
 delete [] used;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new MoleculeCMModifier(args); }
void destroy(Plugin * m) { delete m; }


