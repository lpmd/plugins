//
//
//

#include "rvcorr.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>

using namespace lpmd;

RVCorr::RVCorr(std::string args): Module("rvcorr")
{
 m = NULL;
 AssignParameter("average", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = GetDouble("rcut");
 nb = GetInteger("bins");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 outputfile = GetString("output");
 do_average = GetBool("average");
}

RVCorr::~RVCorr() { if (m != NULL) delete m; }

void RVCorr::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = rvcorr                                                   \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Especifica el numero de divisiones entre 0 y rcut.       \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo de rvcorr     \n";
 std::cout << "      output        : Fichero en el que se graba rvcorr.                       \n";
 std::cout << "      average       : Setea si calculo o no el promedio de cada calculo.       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use rvcorr                                                                    \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     output rvcorr.dat                                                         \n";
 std::cout << "     rcut 15.0                                                                 \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property rvcorr start=1 each=10 end=100                                       \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string RVCorr::Keywords() const { return "rcut bins start end each output average"; }

void RVCorr::Evaluate(SimulationCell & simcell, Potential & pot)
{
 // fabs(rcut) < 1e-05 used to avoid comparing doubles
 if (nb == 0 || fabs(rcut) < 1e-05) throw PluginError("rvcorr", "Error in calculation: Cutoff or bins have wrong value.");
 double dr = rcut/ double(nb);
 int nsp = simcell.SpeciesList().size();
 unsigned long int N = simcell.Size();
 double **g;
 long **nr;
 g = new double*[nb];
 nr = new long*[nb]; 
 for (int i=0;i<nb;i++) 
 { 
  g[i]=new double[(int)(nsp*(nsp+1)/2)]; 
  nr[i]=new long[(int)(nsp*(nsp+1)/2)]; 
 }
 for (int i=0;i<nb;i++) 
 { 
  for (int j=0;j<(int)(nsp*(nsp+1)/2);j++) g[i][j]=0.0e0;
  for (int j=0;j<(int)(nsp*(nsp+1)/2);j++) nr[i][j]=0;
 }
 int s=0;
 const std::list<std::string> lst = simcell.SpeciesPairs();

 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)	   
 {
  //Hace funcional cada una de las especies de los pares.
  std::vector<std::string> loa = SplitTextLine(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  //Cuenta los atomos de cada especie atomica.
  int ne1=0,ne2=0;
  for(unsigned long int m=0;m<N;m++)
  {
   if(simcell[m].Species()==e1) ne1++;
   if(simcell[m].Species()==e2) ne2++;
  }
  //Comienza la iteracion principal
  for (unsigned long int i=0;i<N;++i)
  {
   if (simcell[i].Species()==e1)
   {
    std::vector<Neighbor> nlist;
    simcell.BuildNeighborList (i,nlist,true, simcell.CMCutoff());
    for (unsigned long int j=0;i<nlist.size();++j)
    {
     const Neighbor & nn = nlist[j];
     if(nn.j->Species() == e2)
     {
      if (nn.r <= rcut)
      {
       int ig=(long)floor(nn.r/dr);
       if (fabs(nn.i->Velocity().Mod()*nn.j->Velocity().Mod()) >= FP_ZERO)
       {
        g[ig][s]=g[ig][s]+(Dot(nn.i->Velocity(), nn.j->Velocity()))/(nn.i->Velocity().Mod()*nn.j->Velocity().Mod());
        nr[ig][s]=nr[ig][s]+1;
       }
      }
     }
    }
   }
  }
  s++;
 }
 //
 // Output
 //
 if (m != NULL) delete m;
 m = new Matrix(1 + nsp*(nsp+1)/2, nb);

 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "r");
 int j=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, (*it)+" rvcorr(r)");
  j++;
 }
 // 
 for(int i=0;i<nb;i++)
 {
  m->Set(0, i, dr*i);
  for (int j=0;j<(int)(nsp*(nsp+1)/2);j++) 
  {
   double rv;
   if (nr[i][j] > 0) rv = g[i][j]/nr[i][j];
   else rv = g[i][j];
   m->Set(j+1, i, rv);
  }
 }
 for (int i=0;i<nb;i++) { delete [] g[i]; delete [] nr[i]; }
 delete [] g;
 delete [] nr;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new RVCorr(args); }
void destroy(Module * m) { delete m; }

