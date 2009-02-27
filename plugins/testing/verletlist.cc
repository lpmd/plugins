//
//
//

#include "verletlist.h"

#include <lpmd/simulationcell.h>

#include <cmath>

using namespace lpmd;

VerletListCellManager::VerletListCellManager(std::string args): Module("verletlist")
{ 
 AssignParameter("cutoff", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 cutoff = GetDouble("cutoff");
 oldposition.clear();
 evaluate = true;
 calls=0;
 neigh=NULL;
}

VerletListCellManager::~VerletListCellManager() { if(neigh!=NULL) delete[] neigh; }

void VerletListCellManager::Show(std::ostream & os) const
{
 Module::Show(os);
 if (fabs(cutoff) < 1E-10) os << "   No cutoff was defined." << '\n';
}

std::string VerletListCellManager::Keywords() const { return "cutoff"; }

void VerletListCellManager::Reset() { }

void VerletListCellManager::UpdateCell(SimulationCell & sc) 
{ 
 if (fabs(cutoff) < 1E-10) 
 {
  for (int q=0;q<3;++q)
   if (0.5*sc.GetVector(q).Mod() > cutoff) cutoff = 0.5*sc.GetVector(q).Mod();
 }
}

double VerletListCellManager::Cutoff() const { return cutoff; }

void VerletListCellManager::BuildList(SimulationCell & sc, bool full, double rcu)
{
 if (cutoff == 0) {cutoff = rcu;}
 std::cerr << "cutoff = "  << cutoff << '\n';
 std::cerr << "rcu = " << rcu << '\n';
 if (rcu > cutoff) throw PluginError("verletlist", "Error in cutoff at verletlist, please check that the cutoff is major in plugin");
 const unsigned long n = sc.size();
 if (calls == 0)
 {
  oldposition.reserve(n);
  for(unsigned long int i=0;i<n;++i) oldposition[i]=sc[i].Position();
 }
 if(sc.MetaData().GetInteger("step") % 10 == 0 || evaluate==true)
 {
  std::cerr << "remake the build in step = " << sc.MetaData().GetInteger("step") << '\n';
  evaluate=false;
  for(unsigned long i=0;i<n;++i)
  {
   sc[i].CleanNeighbors();
   if (full)
   {
    //
    for (unsigned long j=0;j<n;++j)
    {
     if (i != j)
     {
      Neighbor nn;
      nn.i = &sc[i];
      nn.j = &sc[j];
      nn.rij = sc.VectorDistance(i, j);
      nn.r = nn.rij.Mod();
      if (cutoff < 1E-30) sc[i].AddNeighbor(nn);
      else if (nn.r < cutoff) sc[i].AddNeighbor(nn);
     }
    }
   }
   else
   {
    //
    for (unsigned long j=i+1;j<n;++j)
    {
     Neighbor nn;
     nn.i = &sc[i];
     nn.j = &sc[j];
     nn.rij = sc.VectorDistance(i, j);
     nn.r = nn.rij.Mod();
     if (cutoff < 1E-30) sc[i].AddNeighbor(nn);
     else if (nn.r < cutoff) sc[i].AddNeighbor(nn);
    }
   }
  }
  //After check if a displacement is more that srcut-frcut.
  double tolerance = cutoff - rcu;
  if(calls!=0)
  {
   for (unsigned long i=0;i<n;i++)
   {
    lpmd::Vector actual = sc[i].Position();
    double displace = (actual - oldposition[i]).Mod();
//    std::cerr << "displace = " << displace << " tolerance = " << tolerance << " -- \n";
    if (displace > tolerance) evaluate=true;
    oldposition[i] = actual;
   }
  }
  calls++;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VerletListCellManager(args); }
void destroy(Module * m) { delete m; }
