/*
 *
 *
 *
 */

#include <lpmd/util.h>
#include "plugincommon.h"
#include "config.h"
#include "version.h"

using namespace lpmd;

const char * PluginVersion()
{
 std::string lver;
 lver = VERSION;
 #ifndef NUMBERED_RELEASE
 lver += " (from ";
 lver += SVNBRANCH;
 lver += (", revision "+ToString<int>(SVNREVISION)+")");
 #endif
 return lver.c_str();
}

lpmd::Matrix* gdr(SimulationCell & simcell,Potential & pot,long int nb,double rcut)
{
 // fabs(rcut) < 1e-05 used to avoid comparing doubles
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
    simcell.BuildNeighborList (i,nlist,true, rcut);
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
 Matrix *m=NULL;
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
 return m;
}
