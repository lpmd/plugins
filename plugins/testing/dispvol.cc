//
//
//

#include "dispvol.h"
#include <lpmd/simulation.h>
#include <lpmd/configurationset.h>

#include <lpmd/properties.h>
#include <lpmd/matrix.h>
#include <lpmd/configuration.h>
#include <lpmd/potential.h>
#include <lpmd/util.h>
#include <lpmd/configurationset.h>
#include <lpmd/storedconfiguration.h>


#include <lpmd/simulation.h>
#include <lpmd/properties.h>

#include <sstream>


#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>


using namespace lpmd;

DispVol::DispVol(std::string args): Module("dispvol")
{
 ParamList & params = (*this);
 AssignParameter("t", "5");
 ProcessArguments(args);
 delta_t = int(params["t"]);
#warning Por qué no pasa outputfile si se inlcuye property.h?
// OutputFile() = params["output"];
}

DispVol::~DispVol() { if (m != NULL) delete m; }

std::string DispVol::Keywords() const { return "t output"; }

void DispVol::Evaluate(lpmd::ConfigurationSet & hist, Potential & pot)
{
 int N = hist.Size();
 int nat = (hist[0].Atoms()).Size();

 //
 // Undo periodicity 
 //
 lpmd::ParticleSet part = hist[0].Atoms();
 lpmd::BasicCell & cell = hist[0].Cell();

 Vector ** noperiodic = new Vector*[N];
 for (int t=0;t<N;++t) noperiodic[t] = new Vector[nat];
 for (int i=0;i<nat;++i) noperiodic[0][i] = part[i].Position();

 for (int t=1;t<N;++t)
  for (int i=0;i<nat;++i)
  {
   part[0].Position() = hist[t-1].Atoms()[i].Position();
   part[1].Position() = hist[t].Atoms()[i].Position();
   noperiodic[t][i] = noperiodic[t-1][i] + cell.Displacement(part[0].Position(), part[1].Position());
  }

 //
 //
 //

 long int t = delta_t;

 lpmd::Matrix & m = CurrentValue();
 m = lpmd::Matrix(2, nat*(N-t));
 m.SetLabel(0, "atom");
 m.SetLabel(1, "dispvol");

 for (int t0=0;t0<N-t;t0++) // loop sobre todos los origenes
 {
  for (int i=0;i<nat;i++)  // loop sobre todos los atomos
  {
   m.Set(0, t0*nat+i, i);
   m.Set(1, t0*nat+i, (M_PI*4.0/3.0)*pow((noperiodic[t0+t][i]-noperiodic[t0][i]).Module(), 3.0));
  }
 }
 for (int i=0;i<N;++i) delete [] noperiodic[i];
 delete [] noperiodic;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new DispVol(args); }
void destroy(Module * m) { delete m; }

