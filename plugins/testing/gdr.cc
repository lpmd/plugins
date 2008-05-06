//
//
//

#include "gdr.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>

using namespace lpmd;

Gdr::Gdr(std::string args): Module("gdr")
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

Gdr::~Gdr() { if (m != NULL) delete m; }

void Gdr::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = gdr                                                      \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular la funcion de autocorrelacion de pares    \n";
 std::cout << " utiliza las condiciones de borde periodicas de la celda.                      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Especifica el numero de divisiones entre 0 y rcut.       \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo de gdr.       \n";
 std::cout << "      output        : Fichero en el que se graba el RDF.                       \n";
 std::cout << "      average       : Setea si calculo o no el promedio de cada calculo.       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use gdr                                                                       \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     output filegdr.dat                                                        \n";
 std::cout << "     rcut 15.0                                                                 \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property gdr start=1 each=10 end=100                                        \n\n";
 std::cout << "      De esta forma clculamos la funcion de distribucion radial de pares en    \n";
 std::cout << " la simulacion entre 1 y 100 cada 10 pasos.                                    \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Gdr::Keywords() const { return "rcut bins start end each output average"; }

void Gdr::Evaluate(SimulationCell & simcell, Potential & pot)
{
 // fabs(rcut) < 1e-05 used to avoid comparing doubles
 if (nb == 0 || fabs(rcut) < 1e-05) throw PluginError("gdr", "Error in calculation: Cutoff or bins have wrong value.");
 double dr = rcut/ double(nb);

 int nsp = simcell.NEspec();
 int N = simcell.Size();
 double **g, *gt;
 g = new double*[nb];
 for(int i=0;i<nb;i++) { g[i]=new double[(int)(nsp*(nsp+1)/2)]; }
 gt = new double[nb]; //total gdr
 for (int i=0;i<nb;i++) 
 { 
  gt[i]=0.0e0;
  for (int j=0;j<(int)(nsp*(nsp+1)/2);j++) g[i][j]=0.0e0;
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
  for(int m=0;m<N;m++)
  {
   if(simcell.GetAtom(m).Species()==e1) ne1++;
   if(simcell.GetAtom(m).Species()==e2) ne2++;
  }
  //Comienza la iteracion principal para el calculo de g(r).
  for(int i=0;i<N;++i)
  {
   if(simcell[i].Species()==e1)
   {
    std::list<Neighbor> nlist;
    simcell.BuildNeighborList (i,nlist,true);
    for(std::list<Neighbor>::const_iterator it=nlist.begin();it!=nlist.end();++it)
    {
     const Neighbor &nn = *it;
     if(nn.j->Species()==e2)
     {
      if(nn.r*nn.r<=rcut*rcut)
      {
       int ig=(long)floor(nn.r/dr);
       g[ig][s]=g[ig][s]+(simcell.Volume())/(4.0e0*M_PI*nn.r*nn.r*dr*ne1*ne2);
      }
     }
    }
   }
  }
  s++;
 }
 //Calcula el valor de g(r) total.
 int j=0;
 for(std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  //Hace funcional cada una de las especies de los pares.
  std::vector<std::string> loa = SplitTextLine(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  //Cuenta la concentracion atomica de cada especie atomica.
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(simcell.GetAtom(m).Species()==e1) ne1++;
   if(simcell.GetAtom(m).Species()==e2) ne2++;
  }
  double ce1 = (double)ne1/(double)N;
  double ce2 = (double)ne2/(double)N;
  //Comienza la asignacion principal para g(r) total.
  for(int i=0;i<nb;i++)
  {
   if(e1==e2) gt[i] = gt[i]+ce1*ce2*g[i][j];
   else {gt[i]=gt[i]+2*ce1*ce2*g[i][j];}
  }
  j++;
 }
 //
 // Output of g(r)
 //
 if (m != NULL) delete m;
 m = new Matrix(2 + nsp*(nsp+1)/2, nb);

 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "r");
 m->SetLabel(nsp*(nsp+1)/2+1, "total g(r)");
 j=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, (*it)+" g(r)");
  j++;
 }
 // 
 for(int i=0;i<nb;i++)
 {
  m->Set(0, i, dr*i);
  for(int j=0;j<(int)(nsp*(nsp+1)/2);j++)
  {
   m->Set(j+1, i, g[i][j]);
  }
  m->Set(nsp*(nsp+1)/2+1, i, gt[i]);
 }
 delete [] gt;
 for (int i=0;i<nb;i++) delete [] g[i];
 delete [] g;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Gdr(args); }
void destroy(Module * m) { delete m; }


