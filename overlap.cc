//
//
//

#include "overlap.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/atompair.h>
#include <lpmd/configuration.h>

#include <sstream>

using namespace lpmd;

Overlap::Overlap(std::string args): Plugin("overlap", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("rcut");
 DefineKeyword("bins", "10");
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("output");
 ProcessArguments(args);
 rcut = params["rcut"];
 bins = int(params["bins"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
}

void Overlap::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo del overlap   \n";
 std::cout << "      bins          : Especifica el numero de intervalos en distancia          \n";
 std::cout << "      output        : Fichero en el que se graba la salida                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use overlap                                                                   \n";
 std::cout << "     output overlap.dat                                                        \n";
 std::cout << "     bins 50                                                                   \n";
 std::cout << "     rcut 15.0                                                                 \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property overlap start=1 each=10 end=100                                    \n\n";
}

double Foverlap(double r0, double r)
{
 if (r > 2*r0) return 0.0;
 else return (1.0 - (1.0/2.0)*(r/r0) - (1.0/8.0)*pow(r/r0, 2.0)+(1.0/16.0)*pow(r/r0, 3.0));
}

void Overlap::Evaluate(Configuration & conf, Potential & pot)
{
 assert(&pot != 0); //icc 869
 BasicParticleSet & atoms = conf.Atoms();
 const long int n = atoms.Size();
 NeighborList * nlistarray = new NeighborList[n];
 for (long int i=0;i<n;++i)
 {
  const NeighborList & nlist = conf.Neighbors(i, true, 2.0*rcut);
  for (long int k=0;k<nlist.Size();++k) nlistarray[i].Append(nlist[k]);
 } 
 //
 // Computation of overlaps
 //
 Matrix & m = CurrentValue();
 m = Matrix(2, bins);
 m.SetLabel(0, "r");
 m.SetLabel(1, "Overlap(r)");
 for (int b=0;b<bins;++b)
 {
  const double r0 = rcut*float(b)/float(bins);
  m.Set(0, b, r0);
  double ov = 0.0;
  long int cnt = 0;
  for (long int i=0;i<n;++i)
  {
   const NeighborList & nlist = nlistarray[i];
   for (int k=0;k<nlist.Size();++k) 
   { 
    ov += Foverlap(r0, sqrt(nlist[k].r2));
    cnt++;
   }
  }
  ov /= float(cnt);
  m.Set(1, b, ov);
 }
 //
 delete [] nlistarray;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Overlap(args); }
void destroy(Plugin * m) { delete m; }

