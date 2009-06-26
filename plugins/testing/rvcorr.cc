//
//
//

#include "rvcorr.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/atompair.h>
#include <lpmd/configuration.h>

#include <sstream>

using namespace lpmd;

RVCorr::RVCorr(std::string args): Plugin("rvcorr", "1.0")
{
 ParamList & params = (*this);
 DefineKeyword("rcut", "10.0");
 DefineKeyword("bins", "200");
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("output");
 DefineKeyword("average", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = double(params["rcut"]);
 nb = int(params["bins"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 do_average = bool(params["average"]);
}

RVCorr::~RVCorr() { }

void RVCorr::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Plugin utilizado para calcular la correlacion de velocidades en funcion de\n";
 std::cout << "      la distancia radial desde un atomo.                                      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Especifica el numero de divisiones entre 0 y rcut.       \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo de rvcorr     \n";
 std::cout << "      output        : Fichero en el que se graba rvcorr.                       \n";
 std::cout << "      average       : Setea si calculo o no el promedio de cada calculo.       \n";
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
}

void RVCorr::Evaluate(Configuration & conf, Potential & pot)
{
 lpmd::BasicParticleSet & atoms = conf.Atoms();
 lpmd::BasicCell & cell = conf.Cell();

 double dr = rcut/ double(nb);
 lpmd::Array<int> esp = atoms.Elements();
 int nsp = esp.Size();
 int N = atoms.Size();

 double **g, *gt;
 g = new double*[nb];
 for(int i=0;i<nb;i++) { g[i]=new double[(int)(nsp*(nsp+1)/2)]; }
 gt = new double[nb]; //total gdr
 for(int i=0;i<nb;i++) 
 { 
  gt[i]=0.0e0;
  for(int j=0;j<(int)(nsp*(nsp+1)/2);j++) g[i][j]=0.0e0;
 }
 int s=0;
 lpmd::Array<std::string> pairs;
 for(int i=0;i<esp.Size();++i)
  for(int j=i;j<esp.Size();++j)
  {
   std::ostringstream ostr;
   ostr << ElemSym[esp[i]]<< "-" << ElemSym[esp[j]];
   pairs.Append(ostr.str());
  }

 for(int i=0;i<pairs.Size();++i)	   
 {
  lpmd::Array<std::string> loa = lpmd::SplitSpeciesPair(pairs[i]); // lista de atomos
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(atoms[m].Symbol()==loa[0]) ne1++;
   if(atoms[m].Symbol()==loa[1]) ne2++;
  }

  for(int i=0;i<N;++i)
  {
   if(atoms[i].Symbol()==loa[0])
   {
    lpmd::NeighborList & nlist = conf.Neighbors(i,true,rcut);
    for(long int k=0; k<nlist.Size();++k)
    {
     const lpmd::AtomPair & nn = nlist[k];
     if(nn.j->Symbol()==loa[1])
     {
      if(nn.r*nn.r<=rcut*rcut)
      {
       int ig=(long)floor(nn.r/dr);
       double rvcorr = (Dot(nn.i->Velocity(), nn.j->Velocity())/(nn.i->Velocity().Module()*nn.j->Velocity().Module()));
       g[ig][s] += rvcorr*(cell.Volume())/(4.0e0*M_PI*nn.r*nn.r*dr*ne1*ne2);
      }
     }
    }
   }
  }
  s++;
 }

 //Calcula el valor de rvcorr(r) total.
 int j=0;
 for(long int i=0;i<pairs.Size();++i)
 {
  lpmd::Array<std::string> loa = lpmd::SplitSpeciesPair(pairs[i]); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(atoms[m].Symbol()==loa[0]) ne1++;
   if(atoms[m].Symbol()==loa[1]) ne2++;
  }
  double ce1 = (double)ne1/(double)N;
  double ce2 = (double)ne2/(double)N;
  for(int i=0;i<nb;i++)
  {
   if(e1==e2) gt[i] = gt[i]+ce1*ce2*g[i][j];
   else {gt[i]=gt[i]+2*ce1*ce2*g[i][j];}
  }
  j++;
 }

 Matrix & m = CurrentValue();
 m = lpmd::Matrix(2 + nsp*(nsp+1)/2, nb);
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "r");
 m.SetLabel(nsp*(nsp+1)/2+1, "total rvcorr(r)");
 j=1;
 for (long int i=0;i<pairs.Size();++i)
 {
  m.SetLabel(j, pairs[i]+" rvcorr(r)");
  j++;
 }
 // 
 for(int i=0;i<nb;i++)
 {
  m.Set(0, i, dr*i);
  for(int j=0;j<(int)(nsp*(nsp+1)/2);j++)
  {
   m.Set(j+1, i, g[i][j]);
  }
  m.Set(nsp*(nsp+1)/2+1, i, gt[i]);
 }
 delete [] gt;
 for (int i=0;i<nb;i++) delete [] g[i];
 delete [] g;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new RVCorr(args); }
void destroy(Plugin * m) { delete m; }

