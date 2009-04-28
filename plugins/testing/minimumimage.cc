//
//
//

#include "minimumimage.h"

#include <lpmd/simulationcell.h>

#include <cmath>

using namespace lpmd;

MinimumImageCellManager::MinimumImageCellManager(std::string args): Module("minimumimage")
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("cutoff", "0.0");
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

void MinimumImageCellManager::Reset() { }

void MinimumImageCellManager::UpdateCell(SimulationCell & sc) 
{ 
 if (fabs(rcut) < 1E-10) 
 {
  for (int q=0;q<3;++q)
   if (0.5*sc.GetVector(q).Module() > rcut) rcut = 0.5*sc.GetVector(q).Module();
 }
}

double MinimumImageCellManager::Cutoff() const { return rcut; }

void MinimumImageCellManager::BuildNeighborList(SimulationCell & sc, long i, std::vector<Neighbor> & nlist, bool full, double rcu)
{
 const long int n = sc.size();
 nlist.clear();
 if (full)
 {
  //
  for (long int j=0;j<n;++j)
  {
   if (i != j)
   {
    Neighbor nn;
    nn.i = &sc[i];
    nn.j = &sc[j];
    nn.rij = sc.VectorDistance(i, j);
    nn.r = nn.rij.Module();
    if (rcut < 1E-30) nlist.push_back(nn);
    else if (nn.r < rcut) nlist.push_back(nn);
   }
  }
 }
 else
 {
  //
  for (long int j=i+1;j<n;++j)
  {
   Neighbor nn;
   nn.i = &sc[i];
   nn.j = &sc[j];
   nn.rij = sc.VectorDistance(i, j);
   nn.r = nn.rij.Module();
   if (rcut < 1E-30) nlist.push_back(nn);
   else if (nn.r < rcut) nlist.push_back(nn);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MinimumImageCellManager(args); }
void destroy(Module * m) { delete m; }


