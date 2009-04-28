//
//
//

#include "localpressure.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>
#include <lpmd/pairpotential.h>
#include <lpmd/potentialarray.h>
#include <lpmd/session.h>

#include <sstream>

using namespace lpmd;

LocalPressure::LocalPressure(std::string args): Module("localpressure")
{
 m = NULL;
 // hasta aqui los valores por omision
 AssignParameter("nx", "10");
 AssignParameter("ny", "10");
 AssignParameter("nz", "10");
 AssignParameter("average", "false");
 ProcessArguments(args);
 n[0] = GetInteger("nx");
 n[1] = GetInteger("ny");
 n[2] = GetInteger("nz");
 rcut = GetDouble("rcut");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 outputfile = GetString("output");
 do_average = GetBool("average");
}

LocalPressure::~LocalPressure() { if (m != NULL) delete m; }

void LocalPressure::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = localpressure                                            \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property localpressure start=1 each=10 end=100                              \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string LocalPressure::Keywords() const { return "nx ny nz rcut start end each output average"; }

void LocalPressure::Evaluate(SimulationCell & simcell, Potential & pot)
{ 
 PotentialArray & p_array = dynamic_cast<PotentialArray &>(pot);

 double cv = simcell.Volume()/(n[0]*n[1]*n[2]);
 if (m != NULL) delete m;
 m = new Matrix(12, n[0]*n[1]*n[2]);
 m->SetLabel(0, "i");
 m->SetLabel(1, "j");
 m->SetLabel(2, "k");
 for (int p=0;p<3;++p)
   for (int q=0;q<3;++q) m->SetLabel(3+(3*p+q), "S"+ToString<int>(p+1)+ToString<int>(q+1));

 int l=0; 
 for (int k=0;k<n[2];++k)
   for (int j=0;j<n[1];++j)
    for (int i=0;i<n[0];++i)
    {
     m->Set(0, l, i);
     m->Set(1, l, j);
     m->Set(2, l, k);
     ++l;
    }

 for (unsigned long int i=0;i<simcell.size();++i)
 {
  int s1, ind[3];
  double stress[3][3];
  for (int p=0;p<3;++p)
    for (int q=0;q<3;++q) stress[p][q]=0.0e0;
  s1 = simcell[i].Species();
  std::vector<Neighbor> nlist;
  simcell.BuildNeighborList(i, nlist, false, rcut);
  Vector fpos = simcell.FracPosition(i);
  for (int q=0;q<3;++q) 
  {
   ind[q] = int(floor(fpos[q]*n[q]));
   if (ind[q] == n[q]) ind[q]--;
  }
  for (std::vector<Neighbor>::const_iterator it=nlist.begin();it!=nlist.end();++it)
  {
   try
   {
    const Neighbor &nn = *it;
    int s2 = (nn.j)->Species();
    PairPotential & ppot = dynamic_cast<PairPotential &>(p_array.Get(s1, s2));
    Vector ff = ppot.pairForce(nn.rij);
    for (int p=0;p<3;++p)
     for (int q=0;q<3;++q) stress[p][q] -= (nn.rij)[p]*ff[q];
   }
   catch (std::exception &e) { throw PluginError("localpressure", "Cannot calculate local stress with a non-pair potential."); }
  }
  int k = ind[0]+n[0]*ind[1]+n[0]*n[1]*ind[2];
  const double pressfactor = GlobalSession.GetDouble("pressfactor");
  for (int p=0;p<3;++p)
   for (int q=0;q<3;++q)
   {
    m->Set(3+(3*p+q), k, m->Get(3+(3*p+q), k)+(pressfactor/cv)*stress[p][q]);
   }
 } 
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new LocalPressure(args); }
void destroy(Module * m) { delete m; }


