//
//
//

#include "localpressure.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <lpmd/pairpotential.h>
#include <lpmd/session.h>
#include <lpmd/basiccell.h>
#include <lpmd/combinedpotential.h>
#include <lpmd/error.h>

#include <sstream>

using namespace lpmd;

LocalPressure::LocalPressure(std::string args): Plugin("localpressure", "2.0")
{
 ParamList & param = (*this);
 // hasta aqui los valores por omision
 AssignParameter("nx", "10");
 AssignParameter("ny", "10");
 AssignParameter("nz", "10");
 ProcessArguments(args);
 n[0] = int(param["nx"]);
 n[1] = int(param["ny"]);
 n[2] = int(param["nz"]);
 rcut = double(param["rcut"]);
 start = int(param["start"]);
 end = int(param["end"]);
 each = int(param["each"]);
 OutputFile() = param["output"];
}

LocalPressure::~LocalPressure() { }

void LocalPressure::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property localpressure start=1 each=10 end=100                              \n\n";
}

void LocalPressure::Evaluate(Configuration & sim, Potential & pot)
{ 
 CombinedPotential * p_array_ptr = 0;
 try { p_array_ptr = dynamic_cast<CombinedPotential *>(&pot); }
 catch (std::exception & ) { throw PluginError("localpressure", "Potentials are not active, cannot compute individual stress tensors"); }
 CombinedPotential & p_array = *p_array_ptr;
 lpmd::BasicParticleSet & atoms = sim.Atoms();
 lpmd::BasicCell & cell = sim.Cell();

 double cv = cell.Volume()/(n[0]*n[1]*n[2]);
 lpmd::Matrix & m = CurrentValue();
 m = lpmd::Matrix(12, n[0]*n[1]*n[2]);
 m.SetLabel(0, "i");
 m.SetLabel(1, "j");
 m.SetLabel(2, "k");
 for (int p=0;p<3;++p)
   for (int q=0;q<3;++q) m.SetLabel(3+(3*p+q), "S"+ToString<int>(p+1)+ToString<int>(q+1));

 int l=0; 
 for (int k=0;k<n[2];++k)
   for (int j=0;j<n[1];++j)
    for (int i=0;i<n[0];++i)
    {
     m.Set(0, l, i);
     m.Set(1, l, j);
     m.Set(2, l, k);
     ++l;
    }

 for (long int i=0;i<atoms.Size();++i)
 {
  int s1, ind[3];
  double stress[3][3];
  for (int p=0;p<3;++p)
    for (int q=0;q<3;++q) stress[p][q]=0.0e0;
  s1 = atoms[i].Z();

  lpmd::NeighborList & nlist = sim.Neighbors(i,false,rcut);
  //simcell.BuildNeighborList(i, nlist, false, rcut);
  Vector fpos = cell.Fractional(atoms[i].Position());
  for (int q=0;q<3;++q) 
  {
   ind[q] = int(floor(fpos[q]*n[q]));
   if (ind[q] == n[q]) ind[q]--;
  }
  for (int j=0;j<nlist.Size();++j)
  {
   try
   {
    const lpmd::AtomPair & nn = nlist[j];
    int s2 = (nn.j)->Z();
    const PairPotential & ppot = dynamic_cast<const PairPotential &>(p_array.PotentialForElements(s1, s2));
    Vector ff = ppot.pairForce(nn.rij);
    for (int p=0;p<3;++p)
     for (int q=0;q<3;++q) stress[p][q] -= (nn.rij)[p]*ff[q];
   }
   catch (std::exception &) { throw PluginError("localpressure", "Cannot calculate local stress with a non-pair potential."); }
  }
  int k = ind[0]+n[0]*ind[1]+n[0]*n[1]*ind[2];
  if (k>=0 && k <= n[0]*n[1]*n[2])
  {
   const double pressfactor = double(GlobalSession["pressfactor"]);
   for (int p=0;p<3;++p)
    for (int q=0;q<3;++q)
    {
     m.Set(3+(3*p+q), k, m.Get(3+(3*p+q), k)+(pressfactor/cv)*stress[p][q]);
    }
  }
 } 
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LocalPressure(args); }
void destroy(Plugin * m) { delete m; }


