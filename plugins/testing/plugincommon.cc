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

lpmd::Matrix* vacf(const std::vector<SimulationCell> & simcell, Potential & Pot, double dt)
{
 int N = simcell.size();
 int nsp = simcell[0].NEspec();
 int *sp = simcell[0].Espec();
 int nat = simcell[0].Size();

 double **vaf=new double*[(int)(N-1)/2];
 for(int i=0;i<(int)(N-1)/2;i++) {vaf[i]=new double[nsp];for(int j=0;j<nsp;j++) vaf[i][j]=0.0e0;}

 Vector ** velocities = new Vector*[N];
 for (int t=0;t<N;++t)  velocities[t] = new Vector[nat];
 if(simcell[0].MetaData().GetInteger("level")==0)
 {
  //
  // Undo periodicity 
  //
  SimulationCell scratch(simcell[0]);
  Vector ** noperiodic = new Vector*[N];
  for (int t=0;t<N;++t) noperiodic[t] = new Vector[nat];
  for (int i=0;i<nat;++i) noperiodic[0][i] = simcell[0].GetAtom(i).Position();
  
  for (int t=1;t<N;++t)
   for (int i=0;i<nat;++i)
   {
    scratch.SetPosition(0, simcell[t-1].GetAtom(i).Position());
    scratch.SetPosition(1, simcell[t].GetAtom(i).Position());
    noperiodic[t][i] = noperiodic[t-1][i] + scratch.VectorDistance(0, 1);
   }
  //
  //Evaluate and set velocities
  //
  for (int i=0;i<nat;++i) velocities[0][i]=(noperiodic[0][i]-noperiodic[N-1][i])/dt;

  for (int t=1;t<N;++t)
   for (int i=0;i<nat;++i)
   {
    Vector vel = (noperiodic[t][i]-noperiodic[t-1][i])/dt;
    velocities[t][i] = vel;
   }
 }
 if(simcell[0].MetaData().GetInteger("level")>=1)
 {
  for (int t=0;t<N;++t)
   for (int i=0;i<nat;++i)
   {
    velocities[t][i] = simcell[t].GetAtom(i).Velocity();
   }
 }
 
 int s=0;
 for(int e1=0;e1<nsp;e1++)	   
 {		 	
   int ne=0;
   for(int i=0;i<nat;i++) {if(simcell[0].GetAtom(i).Species() == sp[e1]) ne++;}
   for(int t0=0;t0<(int)(N-1)/2;t0++)
   {
     for(int t=0;t<(int)(N-1)/2;t++)
     {
      for(int i=0;i<nat;i++)
      {
       Vector v0n = velocities[t0][i];//simcell[t0].GetAtom(i).Velocity();
       Vector v1n = velocities[t0+t][i];//simcell[t0+t].GetAtom(i).Velocity();
       if(simcell[t0].GetAtom(i).Species() == sp[e1])
       {
	  vaf[t][e1]+=Dot(v0n,v1n)/(ne*(int)(N-1)/2);
       }
      }
     }
   }
   s++;
 }

 for (int i=0;i<N;++i) delete[] velocities[i];
 delete [] velocities;

 //
 // Output of vacf
 //
 Matrix *m=NULL;
 m = new Matrix(nsp+1, (int)(N-1)/2);
 const std::list<std::string> lst = simcell[0].SpeciesList();

 int k=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(k, (*it));
  k++;
 }
 m->SetLabel(0,"time");

 for(int i=0;i<(int)(N-1)/2;++i)
 {
  m->Set(0, i, dt*i);
  for (int j=0;j<nsp;++j)
  {
   m->Set(j+1, i, vaf[i][j]);
  }
 }
 for(int i=0;i<(int)(N-1)/2;++i) delete [] vaf[i];
 delete [] vaf;
 return m;
}
