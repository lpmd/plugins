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
 if (rcut <= 0.0) os << "   No cutoff was defined." << '\n';
}

std::string MinimumImageCellManager::Keywords() const { return "cutoff"; }

void MinimumImageCellManager::Reset() { }

void MinimumImageCellManager::UpdateCell(SimulationCell & sc) { }

void MinimumImageCellManager::BuildNeighborList(SimulationCell & sc, long i, std::list<Neighbor> & nlist, bool full)
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
    else if (nn.r < rcut) nlist.push_back(nn);
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
   else if (nn.r < rcut) nlist.push_back(nn);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MinimumImageCellManager(args); }
void destroy(Module * m) { delete m; }


