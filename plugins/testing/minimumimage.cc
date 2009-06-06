//
//
//

#include "minimumimage.h"
#include <lpmd/configuration.h>

#include <cmath>

using namespace lpmd;

MinimumImageCellManager::MinimumImageCellManager(std::string args): Plugin("minimumimage", "2.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("cutoff", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = params["cutoff"];
}

MinimumImageCellManager::~MinimumImageCellManager() { }

void MinimumImageCellManager::Show(std::ostream & os) const
{
 Module::Show(os);
 if (fabs(rcut) < 1E-10) os << "   No cutoff was defined." << '\n';
}

void MinimumImageCellManager::Reset() { }

void MinimumImageCellManager::UpdateCell(Configuration & conf) 
{ 
 BasicCell & cell = conf.Cell();
 if (fabs(rcut) < 1E-10) 
 {
  for (int q=0;q<3;++q)
   if (0.5*cell[q].Module() > rcut) rcut = 0.5*cell[q].Module();
 }
}

double MinimumImageCellManager::Cutoff() const { return rcut; }

void MinimumImageCellManager::BuildNeighborList(Configuration & conf, long i, NeighborList & nlist, bool full, double rcu)
{
 if (rcu < 1.0E-10) EndWithError("MinimumImage cutoff equals zero...");
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 const long int n = atoms.Size();
 nlist.Clear();
 AtomPair nn;
 nn.i = &atoms[i];
 const Vector & ipos = atoms[i].Position(); 
 if (full)
 {
  for (long int j=0;j<n;++j)
   if (i != j)
   {
    nn.j = &(atoms[j]);
    nn.rij = cell.Displacement(ipos, nn.j->Position());
    nn.r = nn.rij.Module();
    if (nn.r < rcu) nlist.Append(nn);
   }
 }
 else
 {
  for (long int j=i+1;j<n;++j)
  {
   nn.j = &(atoms[j]);
   nn.rij = cell.Displacement(ipos, nn.j->Position());
   nn.r = nn.rij.Module();
   if (nn.r < rcu) nlist.Append(nn);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new MinimumImageCellManager(args); }
void destroy(Plugin * m) { delete m; }


