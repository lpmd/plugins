//
//
//

#include "vacf.h"

using namespace lpmd;

Vacf::Vacf(std::string args): Module("vacf")
{
 m = NULL;
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 ProcessArguments(args);
}

Vacf::~Vacf() { if (m != NULL) delete m; }

void Vacf::Evaluate(const std::vector<SimulationCell> & simcell, Potential & pot)
{
 int N = simcell.size();
 std::list<std::string> species = simcell[0].SpeciesList();
 int nsp = species.size(), q = 0;
 int * sp = new int[nsp];
 for (std::list<std::string>::const_iterator it=species.begin();it!=species.end();++it) sp[q++] = ElemNum(*it);
 unsigned long int nat = simcell[0].size();

 double **vaf=new double*[(int)(N-1)/2];
 for(int i=0;i<(int)(N-1)/2;i++) {vaf[i]=new double[nsp];for(int j=0;j<nsp;j++) vaf[i][j]=0.0e0;}

 int s=0;
 for(int e1=0;e1<nsp;e1++)	   
 {		 	
   int ne=0;
   for (unsigned long int i=0;i<nat;i++) { if (simcell[0][i].Species() == sp[e1]) ne++;}
   for(int t0=0;t0<(int)(N-1)/2;t0++)
   {
      for(int t=0;t<(int)(N-1)/2;t++)
      {
       for (unsigned long int i=0;i<nat;i++)
       {
	const Vector & v0n = simcell[t0][i].Velocity();
	Vector v1n = simcell[t0+t][i].Velocity();
	if(simcell[t0][i].Species() == sp[e1]) vaf[t][e1]+=Dot(v0n,v1n)/(ne*(int)(N-1)/2);
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
 delete [] sp;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Vacf(args); }
void destroy(Module * m) { delete m; }
