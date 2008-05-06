//
//
//

#include "pairdistances.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>

using namespace lpmd;

PairDistances::PairDistances(std::string args): Module("pairdistances")
{
 m = NULL;
 ProcessArguments(args);
 rcut = GetDouble("rcut");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 outputfile = GetString("output");
}

PairDistances::~PairDistances() { if (m != NULL) delete m; }

void PairDistances::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = pairdistances                                            \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      rcut          : Especifica el radio maximo para el conteo de pares       \n";
 std::cout << "      output        : Fichero en el que se graba la salida                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use pairdistances                                                             \n";
 std::cout << "     output pd.dat                                                             \n";
 std::cout << "     rcut 15.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property pairdistances start=1 each=10 end=100                                \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string PairDistances::Keywords() const { return "rcut start end step output"; }

void PairDistances::Evaluate(SimulationCell & simcell, Potential & pot)
{
 const long n = simcell.Size();
 std::list<Neighbor> total_list;
 for (long i=0;i<n;++i)
 {
  std::list<Neighbor> nlist;
  simcell.BuildNeighborList(i, nlist, false);
  for(std::list<Neighbor>::const_iterator it=nlist.begin();it!=nlist.end();++it)
  {
   const Neighbor & nn = *it;
   if (nn.r < rcut) total_list.push_back(nn); 
  }
 } 
 //
 // Output 
 //
 if (m != NULL) delete m;
 m = new Matrix(3, total_list.size());
 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "r");
 m->SetLabel(1, "i");
 m->SetLabel(2, "j");
 // 
 long i = 0;
 for (std::list<Neighbor>::const_iterator it=total_list.begin();it!=total_list.end();++it)
 {
  const Neighbor & nn = *it;
  m->Set(0, i, nn.r);
  m->Set(1, i, nn.i->Index());
  m->Set(2, i, nn.j->Index());
  i++;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new PairDistances(args); }
void destroy(Module * m) { delete m; }


