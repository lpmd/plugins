//
//
//

#include "minimumimage.h"

#include <lpmd/simulationcell.h>

#include <cmath>

using namespace lpmd;

MinimumImageCellManager::MinimumImageCellManager(std::string args): Module("minimumimage")
{ 
 AssignParameter("cutoff", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = GetDouble("cutoff");
}

MinimumImageCellManager::~MinimumImageCellManager() { }

void MinimumImageCellManager::Show(std::ostream & os) const
{
 Module::Show(os);
 if (fabs(rcut) < 1E-10) os << "   No cutoff was defined." << '\n';
}

std::string MinimumImageCellManager::Keywords() const { return "cutoff"; }

void MinimumImageCellManager::Reset() { }

void MinimumImageCellManager::UpdateCell(SimulationCell & sc) 
{ 
 if (fabs(rcut) < 1E-10) 
 {
  for (int q=0;q<3;++q)
   if (0.5*sc.GetVector(q).Mod() > rcut) rcut = 0.5*sc.GetVector(q).Mod();
 }
}

double MinimumImageCellManager::Cutoff() const { return rcut; }

void MinimumImageCellManager::BuildNeighborList(SimulationCell & sc, long i, std::list<Neighbor> & nlist, bool full, double rcu)
{
 const long n = sc.Size();
 nlist.clear();
 if (full)
 {
  //
  for (long j=0;j<n;++j)
  {
   if (i != j)
   {
    Neighbor nn;
    nn.i = &sc[i];
    nn.j = &sc[j];
    nn.rij = sc.VectorDistance(i, j);
    nn.r = nn.rij.Mod();
    if (rcut < 1E-30) nlist.push_back(nn);
    else if (nn.r < rcu) nlist.push_back(nn);
   }
  }
 }
 else
 {
  //
  for (long j=i+1;j<n;++j)
  {
   Neighbor nn;
   nn.i = &sc[i];
   nn.j = &sc[j];
   nn.rij = sc.VectorDistance(i, j);
   nn.r = nn.rij.Mod();
   if (rcut < 1E-30) nlist.push_back(nn);
   else if (nn.r < rcu) nlist.push_back(nn);
  }
 }
}

void MinimumImageCellManager::BuildList(SimulationCell &sc, bool full, double rc)
{
 for (long int i=0;i<sc.Size();++i)
 {
  sc.GetAtom(i).CleanNeighbors();
  std::list<Neighbor> tmp;
  BuildNeighborList(sc, i, tmp, full, rc);
  for (std::list<Neighbor>::const_iterator it=tmp.begin();it!=tmp.end();++it) 
      sc.GetAtom(i).Neighbors().Add(*it);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MinimumImageCellManager(args); }
void destroy(Module * m) { delete m; }


