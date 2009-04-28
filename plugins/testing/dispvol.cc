//
//
//

#include "dispvol.h"

using namespace lpmd;

DispVol::DispVol(std::string args): Module("dispvol")
{
 m = NULL;
 AssignParameter("t", "5");
 ProcessArguments(args);
 delta_t = GetInteger("t");
 outputfile = GetString("output");
}

DispVol::~DispVol() { if (m != NULL) delete m; }

std::string DispVol::Keywords() const { return "t output"; }

void DispVol::Evaluate(const std::vector<SimulationCell> & simcell, Potential & pot)
{
 int N = simcell.size();
 int nat = simcell[0].size();

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
 //
 //

 long int t = delta_t;

 m = new Matrix(2, nat*(N-t));
 m->SetLabel(0, "atom");
 m->SetLabel(1, "dispvol");

 for (int t0=0;t0<N-t;t0++) // loop sobre todos los origenes
 {
  for (int i=0;i<nat;i++)  // loop sobre todos los atomos
  {
   m->Set(0, t0*nat+i, i);
   m->Set(1, t0*nat+i, (M_PI*4.0/3.0)*pow((noperiodic[t0+t][i]-noperiodic[t0][i]).Module(), 3.0));
  }
 }
 for (int i=0;i<N;++i) delete [] noperiodic[i];
 delete [] noperiodic;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new DispVol(args); }
void destroy(Module * m) { delete m; }

