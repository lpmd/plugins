//
//
//

#include "pairdistances.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/configuration.h>

#include <sstream>

using namespace lpmd;

PairDistances::PairDistances(std::string args): Module("pairdistances")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("rcut");
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("output");
 ProcessArguments(args);
 rcut = GetDouble("rcut");
 start = GetInteger("start");
 end = GetInteger("end");
 each = GetInteger("each");
 OutputFile() = GetString("output");
}

void PairDistances::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      rcut          : Especifica el radio maximo para el conteo de pares       \n";
 std::cout << "      output        : Fichero en el que se graba la salida                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use pairdistances                                                             \n";
 std::cout << "     output pd.dat                                                             \n";
 std::cout << "     rcut 15.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property pairdistances start=1 each=10 end=100                                \n\n";
}

void PairDistances::Evaluate(Configuration & conf, Potential & pot)
{
 BasicParticleSet & atoms = conf.Atoms();
 const long int n = atoms.Size();
 NeighborList total_list;
 for (long int i=0;i<n;++i)
 {
  const NeighborList & nlist = conf.Neighbors(i, false, rcut);
  for (long int k=0;k<nlist.Size();++k) total_list.Append(nlist[k]);
 } 
 //
 // Output 
 //
 Matrix & m = CurrentValue();
 m = Matrix(3, total_list.Size());
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "r");
 m.SetLabel(1, "i");
 m.SetLabel(2, "j");
 // 
 for (long int i=0;i<total_list.Size();++i)
 {
  const AtomPair & nn = total_list[i];
  m.Set(0, i, nn.r);
  m.Set(1, i, nn.i_index);
  m.Set(2, i, nn.j_index);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new PairDistances(args); }
void destroy(Module * m) { delete m; }


