/*
 *
 *
 *
 */

#include "atomtrail.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/atompair.h>
#include <lpmd/configuration.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace lpmd;

AtomTrail::AtomTrail(std::string args): Plugin("atomtrail", "1.0")
{
 ParamList & params = (*this);

 DefineKeyword("nx", "50");
 DefineKeyword("ny", "50");
 DefineKeyword("plane", "XY");
 DefineKeyword("species", "all");
 
 ProcessArguments(args);

 nx = int(params["nx"]);
 ny = int(params["ny"]);
 plane = params["plane"];
 species = params["species"];
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
}

AtomTrail::~AtomTrail() { }

void AtomTrail::ShowHelp() const
{
 std::cout << " Example :                                                                     \n";
 std::cout << " use atomtrail                                                                 \n";
 std::cout << "     output trail.dat                                                          \n";
 std::cout << "     nx 100                                                                    \n";
 std::cout << "     ny 100                                                                    \n";
 std::cout << "     plane XY                                                                  \n";
 std::cout << " enduse                                                                        \n";
}

void AtomTrail::Evaluate(Configuration & conf, Potential & pot)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 int pl1, pl2;
 if ((plane == "XY") || (plane == "xy")) { pl1 = 0; pl2 = 1; }
 else if ((plane == "XZ") || (plane == "xz")) { pl1 = 0; pl2 = 2; }
 else if ((plane == "YZ") || (plane == "yz")) { pl1 = 1; pl2 = 2; }
 else throw PluginError("atomtrail", "plane must be one of: XY, XZ, YZ");

 double * trail = new double[nx*ny];
 for (int i=0;i<ny;++i)
   for (int j=0;j<nx;++j) trail[j*nx+i] = 0.0;

 long int cnt = 0;
 for (long int i=0;i<atoms.Size();++i)
 {
  if ((species == "all") || (species == atoms[i].Symbol()))
  {
   Vector fpos = cell.Fractional(atoms[i].Position());
   int i = int(floor((nx-1)*fpos[pl1]));
   int j = int(floor((ny-1)*fpos[pl2]));
   trail[j*nx+i] += 1.0;
   cnt++;
  }
 }

 //
 // Output 
 //
 Matrix & m = CurrentValue();
 m = Matrix(3, nx*ny);
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "bin_x");
 m.SetLabel(1, "bin_y");
 m.SetLabel(2, "count");

 long int n=0;
 for (int i=0;i<ny;++i)
   for (int j=0;j<nx;++j)
   {
    m.Set(0, n, j);
    m.Set(1, n, i);
    m.Set(2, n, trail[j*nx+i]/double(cnt));
    n++;
   }

 delete [] trail;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new AtomTrail(args); }
void destroy(Plugin * m) { delete m; }

