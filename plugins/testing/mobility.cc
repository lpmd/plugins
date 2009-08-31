//
//
//

#include "mobility.h"

#include <lpmd/simulationhistory.h>
#include <lpmd/storedconfiguration.h>

using namespace lpmd;

Mobility::Mobility(std::string args): Plugin("mobility", "2.0")
{
 ParamList & params = (*this);
 DefineKeyword("t", "5");
 ProcessArguments(args);
 //
 delta_t = int(params["t"]);
 rcutmin = double(params["rcutmin"]);
 rcutmax = double(params["rcutmax"]);
 outputfile = params["output"];
}

Mobility::~Mobility() { }

void Mobility::Evaluate(ConfigurationSet & hist, Potential & pot)
{
 assert(&pot != 0);//icc 869
 long int N = hist.Size(); // number of configurations, not number of atoms
 DebugStream() << "-> Computing MSD over " << N << " configurations\n";
 long int nat = hist[0].Atoms().Size(); // number of atoms
 DebugStream() << "-> First configuration has " << nat << " atoms\n";

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
   const Vector & v0 = scratch_atoms[0].Position() = hist[t-1].Atoms()[i].Position();
   const Vector & v1 = scratch_atoms[1].Position() = hist[t].Atoms()[i].Position();
   noperiodic[t][i] = noperiodic[t-1][i] + cell.Displacement(v0, v1);
  }

 //
 //
 //
 long int t = delta_t;

 Matrix & m = CurrentValue();
 m = Matrix(5, N-t);
 m.SetLabel(0, "atom");
 m.SetLabel(1, "distance");
 m.SetLabel(2, "J0");
 m.SetLabel(3, "J1");
 m.SetLabel(4, "Jr");

 for (int t0=0;t0<N-t;t0++) // loop sobre todos los origenes
 {
  long int nj0 = 0, nj1 = 0, njr = 0;
  double avdist = 0.0;
  for (int i=0;i<nat;i++)  // loop sobre todos los atomos
  {
   double rr = (noperiodic[t0+t][i]-noperiodic[t0][i]).Module();
   if (rr < rcutmin) nj0++;
   else if (rr > rcutmax) njr++;
   else nj1++;
   avdist += rr;
  }
  avdist /= double(nat);
  m.Set(0, t0, t0);
  m.Set(1, t0, avdist);
  m.Set(2, t0, nj0/double(nj0+nj1+njr));
  m.Set(3, t0, nj1/double(nj0+nj1+njr));
  m.Set(4, t0, njr/double(nj0+nj1+njr));
 }
 for (int i=0;i<N;++i) delete [] noperiodic[i];
 delete [] noperiodic;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Mobility(args); }
void destroy(Plugin * m) { delete m; }

