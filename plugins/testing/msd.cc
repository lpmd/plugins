//
//
//

#include "msd.h"

using namespace lpmd;

MSD::MSD(std::string args): Module("msd")
{
 m = NULL;
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 ProcessArguments(args);
}

MSD::~MSD() { if (m != NULL) delete m; }

void MSD::Evaluate(const std::vector<SimulationCell> & simcell, Potential & pot)
{
 int N = simcell.size();
 std::list<std::string> species = simcell[0].SpeciesList();
 int nsp = species.size(), q = 0;
 int * sp = new int[nsp];
 for (std::list<std::string>::const_iterator it=species.begin();it!=species.end();++it) sp[q++] = ElemNum(*it);
 unsigned long int nat = simcell[0].size();
 double ** msd = new double*[(int)(N-1)/2];
 for (int i=0;i<(int)(N-1)/2;i++) 
 {
  msd[i] = new double[nsp];
  for (int j=0;j<nsp;j++) msd[i][j]=0.0e0;
 }

 //
 // Undo periodicity 
 //
 SimulationCell scratch(simcell[0]);
 Vector ** noperiodic = new Vector*[N];
 for (int t=0;t<N;++t) noperiodic[t] = new Vector[nat];
 for (unsigned long int i=0;i<nat;++i) noperiodic[0][i] = simcell[0][i].Position();
 
 for (int t=1;t<N;++t)
  for (unsigned long int i=0;i<nat;++i)
  {
   scratch.SetPosition(0, simcell[t-1][i].Position());
   scratch.SetPosition(1, simcell[t][i].Position());
   noperiodic[t][i] = noperiodic[t-1][i] + scratch.VectorDistance(0, 1);
  }

 //
 //
 //

 int s=0;
 for(int e1=0;e1<nsp;e1++)
 {		 	
  int ne = 0;
  for (unsigned long int i=0;i<nat;i++) 
  {
   if (simcell[0][i].Species() == sp[e1]) ne++;
  }
  for (int t0=0;t0<(int)(N-1)/2;t0++) // loop sobre todos los origenes
  {
   for (int t=0;t<(int)(N-1)/2;t++) // loop sobre la separacion en tiempo
   {
    for (unsigned long int i=0;i<nat;i++)  // loop sobre todos los atomos
    {
     if (simcell[t0][i].Species() == sp[e1]) msd[t][e1] += (noperiodic[t0+t][i]-noperiodic[t0][i]).SquareModule();
    }
   }
  }
  s++;
 }
 for (int i=0;i<N;++i) delete [] noperiodic[i];
 delete [] noperiodic;

 //
 // Output MSD
 //
 m = new Matrix(nsp+1, (int)(N-1)/2);
 m->SetLabel(0, "time");
 for (int j=0;j<nsp;++j) m->SetLabel(j+1, "MSD");
 for(int i=0;i<(int)(N-1)/2;++i)
 {
  m->Set(0, i, i);
  for(int j=0;j<nsp;++j) m->Set(j+1, i, msd[i][j]/double(nat*(N-1)/2));
 }
 delete [] sp;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MSD(args); }
void destroy(Module * m) { delete m; }

