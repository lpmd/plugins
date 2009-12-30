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
 CombinedPotential & comb = dynamic_cast<CombinedPotential &>(pot);
 for (long i=0;i<n;++i) pair_energy[i] = comb.AtomEnergy(conf, i);
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

