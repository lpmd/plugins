//
//
//


#include <sstream>

#include "vacf.h"

using namespace lpmd;

Vacf::Vacf(std::string args): Module("vacf")
{
 m = NULL;
 ProcessArguments(args);
}

Vacf::Vacf(double rc, int b): Module("vacf")
{
 m = NULL;
}

Vacf::~Vacf()
{
 if (m != NULL) delete m;
}

void Vacf::Evaluate(const std::vector<SimulationCell> & simcell, Potential & pot)
{
 int N = simcell.size();
 int nsp = simcell[0].NEspec();
 int *sp = simcell[0].Espec();
 int nat = simcell[0].Size();

 double **vaf=new double*[(int)(N-1)/2];
 for(int i=0;i<(int)(N-1)/2;i++) {vaf[i]=new double[nsp];for(int j=0;j<nsp;j++) vaf[i][j]=0.0e0;}

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
	Vector v0n = simcell[t0].GetAtom(i).Velocity();
	Vector v1n = simcell[t0+t].GetAtom(i).Velocity();
	if(simcell[t0].GetAtom(i).Species() == sp[e1])
	{
	   vaf[t][e1]+=Dot(v0n,v1n)/(ne*(int)(N-1)/2);
	}
       }
      }
   }
   s++;
 }

 //
 // Output of vacf
 //
 m = new Matrix(nsp, (int)(N-1)/2);

 for(int i=0;i<(int)(N-1)/2;i++)
 {
//  m->Set(0, i, dr*i);
  for(int j=0;j<nsp;j++)
  {
   m->Set(j+1, i, vaf[i][j]);
  }
 }
}

void Vacf::SetParameter(std::string name) { }

void Vacf::Show() const
{
 Module::Show();
}

std::string Vacf::Keywords() const { return ""; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Vacf(args); }
void destroy(Module * m) { delete m; }


