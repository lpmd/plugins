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

//---------------------------------------------------------------//
// AZAR
//
// 11Nov99 primera version
// 30Nov99 header externo y sobrecarga de randomize
//________________________
//JR//

void randomize(unsigned int seed)
{
  srandom(seed);	// Inicializa el random
}

void randomize()
{
  srandom(time(0) * getpid());	// Inicializa el random
}
  
double dazar(double inf, double sup)
{
  return (sup-inf)* double(random())/double(0x7fffffff)+inf ;
}

int iazar(int inf, int sup)
{
  return int(floor(dazar(inf, sup+1))) ;
}
//---------------------------------------------------------------//


lpmd::Matrix* gdr(SimulationCell & simcell,Potential & pot,long int nb,double rcut)
{
 // fabs(rcut) < 1e-05 used to avoid comparing doubles
 double dr = rcut/ double(nb);

 int nsp = simcell.SpeciesList().size();
 int N = simcell.size();
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
  std::vector<std::string> loa = StringSplit< std::vector<std::string> >(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  //Cuenta los atomos de cada especie atomica.
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(simcell[m].Species()==e1) ne1++;
   if(simcell[m].Species()==e2) ne2++;
  }
  //Comienza la iteracion principal para el calculo de g(r).
  for(int i=0;i<N;++i)
  {
   if(simcell[i].Species()==e1)
   {
    std::vector<Neighbor> nlist;
    simcell.BuildNeighborList(i, nlist, true, rcut);
    for(std::vector<Neighbor>::const_iterator it=nlist.begin();it!=nlist.end();++it)
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
  std::vector<std::string> loa = StringSplit< std::vector<std::string> >(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  //Cuenta la concentracion atomica de cada especie atomica.
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(simcell[m].Species()==e1) ne1++;
   if(simcell[m].Species()==e2) ne2++;
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
 int nsp = simcell[0].SpeciesList().size();
 std::list<std::string> species = simcell[0].SpeciesList();
 int * sp = new int[nsp];
 int q = 0;
 for (std::list<std::string>::const_iterator it=species.begin();it!=species.end();++it) sp[q++] = ElemNum(*it);
 int nat = simcell[0].size();

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
  for (int i=0;i<nat;++i) noperiodic[0][i] = simcell[0][i].Position();
  
  for (int t=1;t<N;++t)
   for (int i=0;i<nat;++i)
   {
    scratch.SetPosition(0, simcell[t-1][i].Position());
    scratch.SetPosition(1, simcell[t][i].Position());
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
    velocities[t][i] = simcell[t][i].Velocity();
   }
 }
 
 int s=0;
 for(int e1=0;e1<nsp;e1++)	   
 {		 	
   int ne=0;
   for(int i=0;i<nat;i++) {if(simcell[0][i].Species() == sp[e1]) ne++;}
   for(int t0=0;t0<(int)(N-1)/2;t0++)
   {
     for(int t=0;t<(int)(N-1)/2;t++)
     {
      for(int i=0;i<nat;i++)
      {
       Vector v0n = velocities[t0][i];//simcell[t0][i].Velocity();
       Vector v1n = velocities[t0+t][i];//simcell[t0+t][i].Velocity();
       if(simcell[t0][i].Species() == sp[e1])
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

void Replicate(SimulationCell & sc, unsigned long nx, unsigned long ny, unsigned long nz)
{
 int Ntmp = sc.size();
 Atom *atomos;
 atomos = new Atom[Ntmp];
 for(int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 // FIXME from 0.5: sc.Initialize(nx*Ntmp);
 for(int i=0;i<Ntmp;i++){sc.Create(new Atom(atomos[i]));}

 for(unsigned long i=1;i<nx;i++)
 {
  for(int j=0;j<Ntmp;j++)
  {
   Atom tmp=atomos[j];
   tmp.SetPos(atomos[j].Position()+sc.GetCell()[0]*i);
   sc.Create(new Atom(tmp));
  }
 }
 delete[] atomos;
 Ntmp = sc.size();
 atomos = new Atom[Ntmp];
 for(int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 // FIXME from 0.5: sc.Initialize(ny*Ntmp);
 for(int i=0;i<Ntmp;i++){sc.Create(new Atom(atomos[i]));}
 for(unsigned long i=1;i<ny;i++)
 {
  for(int j=0;j<Ntmp;j++)
  {
   Atom tmp=atomos[j];
   tmp.SetPos(atomos[j].Position()+sc.GetCell()[1]*i);
   sc.Create(new Atom(tmp));
  }
 }
 delete[] atomos;
 Ntmp = sc.size();
 atomos = new Atom[Ntmp];
 for(int i=0;i<Ntmp;i++) atomos[i]=sc[i];
 // FIXME from 0.5: sc.Initialize(nz*Ntmp);
 for(int i=0;i<Ntmp;i++) { sc.Create(new Atom(atomos[i]));}
 for(unsigned long i=1;i<nz;i++)
 {
  for(int j=0;j<Ntmp;j++)
  {
   Atom tmp=atomos[j];
   tmp.SetPos(atomos[j].Position()+sc.GetCell()[2]*i);
   sc.Create(new Atom(tmp));
  }
 }
 delete[] atomos;
 //Resetea los vectores base de la celda.
 Vector a=sc.GetCell()[0];
 sc.GetCell()[0] = a*nx;
 Vector b=sc.GetCell()[1];
 sc.GetCell()[1] = b*ny;
 Vector c=sc.GetCell()[2];
 sc.GetCell()[2] = c*nz;
 //Asigna el index() a cada atomo de la celda.
 sc.AssignIndex();
 sc.ClearForces();
}

void Rotate(SimulationCell & sc, lpmd::Vector rotate)
{
 // euler rotation matrix
 double rotmat[3][3];
 // Eulerian Angles
 double phi=rotate[0], psi=rotate[1], theta=rotate[2];
// double phi=M_PI/8.4, psi=M_PI/8.4, theta=0;
 // Eulerian Matrix
 rotmat[0][0] = (cos(phi)*cos(psi)-cos(theta)*sin(phi)*sin(psi));
 rotmat[0][1] = (cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi));
 rotmat[0][2] = (sin(theta)*sin(psi));
 rotmat[1][0] = (-cos(theta)*cos(psi)*sin(phi)-cos(phi)*sin(psi));
 rotmat[1][1] = ( cos(theta)*cos(phi)*cos(psi)-sin(phi)*sin(psi));
 rotmat[1][2] = (cos(psi)*sin(theta));
 rotmat[2][0] = (sin(theta)*sin(phi));
 rotmat[2][1] = (-cos(phi)*sin(theta));
 rotmat[2][2] = ( cos(theta));

 // We rotate the simulation cell vectors
 for (int n=0; n<3; n++)
 {
  double v1=sc.GetCell()[n][0], v2=sc.GetCell()[n][1], v3=sc.GetCell()[n][2];
  Vector newvec;
  for (int i=0; i<3; i++) newvec[i] = v1*rotmat[i][0]+v2*rotmat[i][1]+v3*rotmat[i][2];
  sc.GetCell()[n] = newvec;
 }
 // We rotate the atoms of the cell
// std::cerr <<"nuevos vectores base de la celda replicada:\n";
// std::cerr << sc.GetVector(0) <<"\n";
// std::cerr << sc.GetVector(1) <<"\n";
// std::cerr << sc.GetVector(2) <<"\n";
 Vector center = (sc.GetCell()[0]+sc.GetCell()[1]+sc.GetCell()[2])*0.5;
// std::cerr << "centro en "<<center <<"\n";
 for (unsigned long n=0;n<sc.size();n++)
 {
//  std::cerr <<"rotando los atomos de la celda replicada.\n";
  Vector pos0 = sc[n].Position();
  double v1=pos0[0], v2=pos0[1], v3=pos0[2];
  Vector newvec;
  for (int i=0; i<3; i++) newvec[i] = v1*rotmat[i][0]+v2*rotmat[i][1]+v3*rotmat[i][2];
//  std::cerr <<"antes de rotar, pos.atom "<<n<<" = "<<pos0<<" de componentes "<<v1<<v2<<v3<<"\n";
  sc.SetPosition(n, newvec);
//  std::cerr <<"despues de rotar, newvec = pos.atom("<<n<<")= "<<newvec<<"\n";
//  std::cerr <<"asignando lo anterior al atom "<<n<<" = "<<sc[n].Position()<<", .\n";
 }
}

void ReplicateRotate(const SimulationCell basecell, lpmd::Vector &cellcenter, lpmd::Vector &CellColor, unsigned long na, unsigned long nb, unsigned long nc, lpmd::Vector rotate, SimulationCell & simcell)
{
 SimulationCell tmpSC=basecell;
 Replicate(tmpSC,na,nb,nc);
/* // CONVERT tmpSC IN A CIRCLE
 for (long i=0; i<tmpSC.size(); i++)
 {
  Vector dist=0.5*(tmpSC.GetVector(0)+tmpSC.GetVector(1)+tmpSC.GetVector(2))-tmpSC[i].Position();
  if (dist.Mod()>0.4*tmpSC.GetVector(0).Mod())
  {
   tmpSC.DeleteAtom(i);  i--;
  }
 }
*/
 Rotate(tmpSC, rotate);
 Vector centro=0.5*(tmpSC.GetCell()[0]+tmpSC.GetCell()[1]+tmpSC.GetCell()[2]);
 unsigned long N = tmpSC.size();
 Atom *atomos = new Atom[N];
 double r=dazar(0,1), g=dazar(0,1), b=dazar(0,1);
 CellColor=r*e1+g*e2+b*e3;
 for(unsigned long i=0;i<N;i++)
 {
  atomos[i]=tmpSC[i];
  Vector vct=atomos[i].Position()+(cellcenter-centro);
  Atom tmp(atomos[i].Species(),vct);
  tmp.SetColor(r*e1+g*e2+b*e3);
  simcell.Create(new Atom(tmp));
 }
 delete[] atomos;
}
