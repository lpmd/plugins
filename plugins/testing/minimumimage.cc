//
//
//

#include "minimumimage.h"

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

void MinimumImageCellManager::UpdateCell(BasicParticleSet & atoms, BasicCell & cell) 
{ 
 if (fabs(rcut) < 1E-10) 
 {
  for (int q=0;q<3;++q)
   if (0.5*cell[q].Module() > rcut) rcut = 0.5*cell[q].Module();
 }
}

double MinimumImageCellManager::Cutoff() const { return rcut; }

void MinimumImageCellManager::BuildNeighborList(BasicParticleSet & atoms, BasicCell & cell, long i, Array<Neighbor> & nlist, bool full, double rcu)
{
 const long int n = atoms.Size();
 nlist.Clear();
 if (full)
 {
  //
  for (long int j=0;j<n;++j)
   if (i != j)
   {
    Neighbor nn;
    nn.i = &atoms[i];
    nn.j = &atoms[j];
    nn.rij = cell.Displacement(atoms[i].Position(), atoms[j].Position());
    nn.r = nn.rij.Module();
    if (rcu < 1E-30) nlist.Append(nn);
    else if (nn.r < rcu) nlist.Append(nn);
   }
 }
 else
 {
  //
  for (long int j=i+1;j<n;++j)
  {
   Neighbor nn;
   nn.i = &atoms[i];
   nn.j = &atoms[j];
   nn.rij = cell.Displacement(atoms[i].Position(), atoms[j].Position());
   nn.r = nn.rij.Module();
   if (rcu < 1E-30) nlist.Append(nn);
   else if (nn.r < rcu) nlist.Append(nn);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MinimumImageCellManager(args); }
void destroy(Module * m) { delete m; }


