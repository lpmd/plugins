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

DispVol::DispVol(std::string args): Plugin("dispvol", "2.0")
{
 ParamList & params = (*this);
 DefineKeyword("t");
 DefineKeyword("mode", "volume");
 DefineKeyword("bins", "300");
 DefineKeyword("maxlength", "15.0");
 ProcessArguments(args);
 delta_t = int(params["t"]);
 mode = params["mode"];
 nbins = int(params["bins"]);
 maxlength = double(params["maxlength"]);
}

DispVol::~DispVol() { if (m != NULL) delete m; }

void DispVol::Evaluate(lpmd::ConfigurationSet & hist, Potential & pot)
{
 assert(&pot != 0);//icc869
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
 {
  DebugStream() << "-> MSD: Undoing periodicity, configuration " << t << '\n';
  for (int i=0;i<nat;++i)
  {
   const Vector v1 = hist[t-1].Atoms()[i].Position();
   const Vector v2 = hist[t].Atoms()[i].Position();
   noperiodic[t][i] = noperiodic[t-1][i] + cell.Displacement(v1, v2);
  }
 }

 //
 //
 //

 long int t = delta_t;

 lpmd::Matrix & m = CurrentValue();
 if (mode == "volume") m = lpmd::Matrix(2, nat*(N-t));
 else if (mode == "vector") m = lpmd::Matrix(4, nat*(N-t));
 else if (mode == "vectorhistogram") m = lpmd::Matrix(4, nbins);
 m.SetLabel(0, "atom");
 if (mode == "volume") m.SetLabel(1, "dispvol");
 else if (mode == "vector")
 {
  m.SetLabel(1, "dx");
  m.SetLabel(2, "dy");
  m.SetLabel(3, "dz");
 }
 else if (mode == "vectorhistogram")
 {
  m.SetLabel(0, "X/Y/Z");
  m.SetLabel(1, "P(X)");
  m.SetLabel(2, "P(Y)");
  m.SetLabel(3, "P(Z)");
 }

 for (int t0=0;t0<N-t;t0++) // loop sobre todos los origenes
 {
  if (mode == "volume")
  {
   for (int i=0;i<nat;i++)  // loop sobre todos los atomos
   {
    m.Set(0, t0*nat+i, i);
    m.Set(1, t0*nat+i, (M_PI*4.0/3.0)*pow((noperiodic[t0+t][i]-noperiodic[t0][i]).Module(), 3.0));
   }
  }
  else if (mode == "vector")
  {
   Vector delta_r;
   for (int i=0;i<nat;i++)  // loop sobre todos los atomos
   {
    m.Set(0, t0*nat+i, i);
    delta_r = noperiodic[t0+t][i]-noperiodic[t0][i];
    for (int q=0;q<3;++q) m.Set(q+1, t0*nat+i, delta_r[q]);
   }
  }
 }
 if (mode == "vectorhistogram")
 {
  double * histx = new double[nbins];
  double * histy = new double[nbins];
  double * histz = new double[nbins];
  int k;
  Vector delta_r;
  for (int i=0;i<nbins;++i) 
      histx[i] = histy[i] = histz[i] = 0.0;
  for (int t0=0;t0<N-t;t0++)
    for (int i=0;i<nat;++i)
    {
     delta_r = noperiodic[t0+t][i]-noperiodic[t0][i];
     k = long(floor(nbins*0.5*((delta_r[0]/maxlength)+1.0)));
     if (k < 0) k = 0;
     else if (k > nbins-1) k = nbins-1;
     histx[k] += 1.0; 
     k = long(floor(nbins*0.5*((delta_r[1]/maxlength)+1.0)));
     if (k < 0) k = 0;
     else if (k > nbins-1) k = nbins-1;
     histy[k] += 1.0; 
     k = long(floor(nbins*0.5*((delta_r[2]/maxlength)+1.0)));
     if (k < 0) k = 0;
     else if (k > nbins-1) k = nbins-1;
     histz[k] += 1.0; 
    }
  for (int i=0;i<nbins;++i)
  {
   m.Set(0, i, (2.0*(double(i)/nbins)-1.0)*maxlength);
   m.Set(1, i, histx[i]/(double(nat*(N-t))));
   m.Set(2, i, histy[i]/(double(nat*(N-t))));
   m.Set(3, i, histz[i]/(double(nat*(N-t))));
  }
  delete [] histx;
  delete [] histy;
  delete [] histz;
 }
 for (int i=0;i<N;++i) delete [] noperiodic[i];
 delete [] noperiodic;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new DispVol(args); }
void destroy(Plugin * m) { delete m; }

