//
//
//

#include "msd.h"
#include <lpmd/simulationhistory.h>
#include <lpmd/storedconfiguration.h>

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

void MSD::Evaluate(const SimulationHistory & hist, Potential & pot)
{
 long int N = hist.Size(); // number of configurations, not number of atoms
 const Array<int> & elements = hist[0].Atoms().Elements();
 int nsp = elements.Size();
 long int nat = hist[0].Atoms().Size(); // number of atoms
 double ** msd = new double*[(int)(N-1)/2];
 for (int i=0;i<(int)(N-1)/2;i++) 
 {
  msd[i] = new double[nsp];
  for (int j=0;j<nsp;j++) msd[i][j]=0.0e0;
 }

 //
 // Undo periodicity 
 //
 StoredConfiguration scratch(hist[0]);
 Vector ** noperiodic = new Vector*[N];
 for (int t=0;t<N;++t) noperiodic[t] = new Vector[nat];
 BasicParticleSet & scratch_atoms = scratch.Atoms();
 BasicCell & cell = scratch.Cell();
 for (long int i=0;i<nat;++i) noperiodic[0][i] = scratch_atoms[i].Position();
 for (int t=1;t<N;++t)
  for (long int i=0;i<nat;++i)
  {
   scratch_atoms[0].Position() = hist[t-1].Atoms()[i].Position();
   scratch_atoms[1].Position() = hist[t].Atoms()[i].Position();
   // FIXME: inlining por mientras tenemos un metodo como sc.VectorDistance
   const Vector & v0 = scratch_atoms[0].Position();
   const Vector & v1 = scratch_atoms[1].Position();
   noperiodic[t][i] = noperiodic[t-1][i] + cell.Displacement(v0, v1);
   // noperiodic[t][i] = noperiodic[t-1][i] + scratch.VectorDistance(0, 1);
  }

 //
 //
 //

 int s=0;
 for(int e1=0;e1<nsp;e1++)
 {		 	
  int ne = 0;
  for (long int i=0;i<nat;i++) 
  {
   if (hist[0].Atoms()[i].Z() == elements[e1]) ne++;
  }
  for (int t0=0;t0<(int)(N-1)/2;t0++) // loop sobre todos los origenes
  {
   for (int t=0;t<(int)(N-1)/2;t++) // loop sobre la separacion en tiempo
   {
    for (long int i=0;i<nat;i++)  // loop sobre todos los atomos
    {
     if (hist[t0].Atoms()[i].Z() == elements[e1]) msd[t][e1] += (noperiodic[t0+t][i]-noperiodic[t0][i]).SquareModule();
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
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new MSD(args); }
void destroy(Module * m) { delete m; }

