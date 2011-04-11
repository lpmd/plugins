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
 DefineKeyword("nz", "50");
 DefineKeyword("mode", "2D");
 DefineKeyword("plane", "XY");
 DefineKeyword("species", "all");
 DefineKeyword("debug", "none");
 
 ProcessArguments(args);

 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);
 mode = params["mode"];
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
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin determine the atom-trail of the atoms in the simulation       \n";
 std::cout << " the atom-trail could be stored in 2D and 3D format.                         \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      e             : The value of epsilon in the potential. [eV]              \n";
 std::cout << "      a             : The value of the a constant in the potencial. [A]        \n";
 std::cout << "      n             : Value of n in the potential.  [integer]                  \n";
 std::cout << "      m             : Value of m in the potential.  [integer]                  \n";
 std::cout << "      c             : Value of the c constant of the potential. [real]         \n";
 std::cout << "      cutoff        : Cutoff of the interatomic potential. [A]                 \n";
 std::cout << "      corrections   : Include(true) or not(false/default) corrections to the   \n";
 std::cout << "                      metallic potential (recommend for homogeneous systems).  \n";
 std::cout << '\n';
 std::cout << " Example          >>                                                           \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use atomtrail as at                                                      \n";
 std::cout << "     output trail2d.dat                                                        \n";
 std::cout << "     nx 100                                                                    \n";
 std::cout << "     ny 100                                                                    \n";
 std::cout << "     plane XY                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Using the loaded plugin :                                                    \n";
 std::cout << " property at start=0 end=-1 each=1                                             \n";
}

void AtomTrail::Evaluate(Configuration & conf, Potential & pot)
{
 assert(&pot != 0);//icc 869
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 int pl1 = -1, pl2 = -1;
 if (mode == "2D")
 {
  if ((plane == "XY") || (plane == "xy")) { pl1 = 0; pl2 = 1; }
  else if ((plane == "XZ") || (plane == "xz")) { pl1 = 0; pl2 = 2; }
  else if ((plane == "YZ") || (plane == "yz")) { pl1 = 1; pl2 = 2; }
  else throw PluginError("atomtrail", "plane must be one of: XY, XZ, YZ");
 }

 double * trail = 0;
 if (mode == "2D")
 { 
  trail = new double[nx*ny]; 
  for (int i=0;i<ny;++i)
    for (int j=0;j<nx;++j) trail[i*nx+j] = 0.0;
 }
 else if (mode == "3D") 
 {
  trail = new double[nx*ny*nz]; 
  for (int i=0;i<nz;++i)
    for (int j=0;j<ny;++j)
      for (int k=0;k<nx;++k) trail[i*ny*nx+j*nx+k] = 0.0;
 }
 else throw PluginError("atomtrail", "Unknown mode, "+mode);

 long int cnt = 0;
 for (long int a=0;a<atoms.Size();++a)
 {
  if ((species == "all") || (species == atoms[a].Symbol()))
  {
   Vector fpos = cell.Fractional(atoms[a].Position());
   if (mode == "2D")
   {
    int i = int(floor((ny-1)*fpos[pl2]));
    int j = int(floor((nx-1)*fpos[pl1]));
    trail[i*nx+j] += 1.0;
   }
   else
   {
    int i = int(floor((nz-1)*fpos[2]));
    int j = int(floor((ny-1)*fpos[1]));
    int k = int(floor((nx-1)*fpos[0]));
    trail[i*ny*nx+j*nx+k] += 1.0;
   }
   cnt++;
  }
 }

 //
 // Output 
 //
 Matrix & m = CurrentValue();
 if (mode == "2D") m = Matrix(3, nx*ny);
 else m = Matrix(4, nx*ny*nz);

 m.SetLabel(0, "bin_x");
 m.SetLabel(1, "bin_y");
 if (mode == "2D") m.SetLabel(2, "count");
 else { m.SetLabel(2, "bin_z"); m.SetLabel(3, "count"); }

 long int n = 0;
 if (mode == "2D")
 {
  for (int i=0;i<ny;++i)
    for (int j=0;j<nx;++j)
    {
     m.Set(0, n, j);
     m.Set(1, n, i);
     m.Set(2, n, trail[i*nx+j]/double(cnt));
     n++;
    }
 }
 else
 {
  for (int i=0;i<nz;++i)
    for (int j=0;j<ny;++j)
      for (int k=0;k<nx;++k)
      {
       m.Set(0, n, k);
       m.Set(1, n, j);
       m.Set(2, n, i);
       m.Set(3, n, trail[i*ny*nx+j*nx+k]/double(cnt));
       n++;
      }
 }
 delete [] trail;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new AtomTrail(args); }
void destroy(Plugin * m) { delete m; }

