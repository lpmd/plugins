//
//
//

#include "pairenergy.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/atompair.h>
#include <lpmd/configuration.h>
#include <lpmd/pairpotential.h>
#include <lpmd/combinedpotential.h>

#include <sstream>
#include <exception>

using namespace lpmd;

PairEnergy::PairEnergy(std::string args): Plugin("pairenergy", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("output");
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
}

void PairEnergy::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular la energia potencial por atomo, para el caso\n"; 
 std::cout << "      de un potencial de pares (por ej. Lennard-Jones).                        \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      output        : Fichero en el que se graba la salida                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use pairenergy                                                                \n";
 std::cout << "     output energ.dat                                                          \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property pairenergy start=1 each=10 end=100                                   \n\n";
}

void PairEnergy::Evaluate(Configuration & conf, Potential & pot)
{
 assert(&pot != 0);
 BasicParticleSet & atoms = conf.Atoms();
 const long int n = atoms.Size();
 double * pair_energy = new double[n];
 for (long i=0;i<n;++i) pair_energy[i] = 0.0;
 const CombinedPotential & comb = dynamic_cast<CombinedPotential &>(pot);
 double rcut = 0.0;
 for (int q=0;q<comb.Size();++q)
     if (comb[q].GetCutoff() > rcut) rcut = comb[q].GetCutoff();
 try
 {
  // Construct an "index table" so we don't have to depend on Atom::Index()
  std::map<BasicAtom *, long int> indices;
  for (long int i=0;i<n;++i) indices[&atoms[i]] = i;

  for (long int i=0;i<n;++i)
  {
   const NeighborList & nlist = conf.Neighbors(i, false, rcut);
   for (long int k=0;k<nlist.Size();++k) 
   {
    const AtomPair & nn = nlist[k];
    int s1 = atoms[i].Z();
    int s2 = nn.j->Z();
    const PairPotential & ppot = dynamic_cast<const PairPotential &>(comb.PotentialForElements(s1, s2));
    const double rcut = ppot.GetCutoff();
    if (nn.r2 >= rcut*rcut) continue;
    else
    {
     const double pe = ppot.pairEnergy(sqrt(nn.r2));
     pair_energy[i] += (0.5*pe);
     pair_energy[indices[nn.j]] += (0.5*pe);
    }
   }
  }
 }
 catch (std::exception & e)
 {
  throw PluginError("pairenergy", "This plugin only makes sense with pair potentials.\n");
 }
 //
 // Output 
 //
 Matrix & m = CurrentValue();
 m = Matrix(2, n);
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "index");
 m.SetLabel(1, "energy");
 // 
 for (long int i=0;i<n;++i)
 {
  m.Set(0, i, i);
  m.Set(1, i, pair_energy[i]);
 }
 delete [] pair_energy;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new PairEnergy(args); }
void destroy(Plugin * m) { delete m; }

