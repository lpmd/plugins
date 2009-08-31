//
//
//

#include "sitecoord.h"

#include <lpmd/simulation.h>
#include <lpmd/properties.h>

#include <sstream>

using namespace lpmd;

SiteCoord::SiteCoord(std::string args): Plugin("sitecoord", "2.0")
{
 ParamList & params=(*this);
 //
 DefineKeyword("rcut");
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("output");
 DefineKeyword("average", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = params["rcut"];
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 do_average = bool(params["average"]);
}

SiteCoord::~SiteCoord() { }

void SiteCoord::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular el numero de coordinacion para atomos     \n";
 std::cout << "      individuales.                                                            \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      rcut          : Especifica el radio maximo para el calculo del numero de coordinacion\n";
 std::cout << "      output        : Fichero en el que se graba la salida.                    \n";
 std::cout << "      average       : Setea si calculo o no el promedio de cada calculo.       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use sitecoord                                                                 \n";
 std::cout << "     output sitecoord.dat                                                      \n";
 std::cout << "     rcut 4.0                                                                 \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property sitecoord start=1 each=10 end=100                                    \n\n";
 std::cout << "      De esta forma calculamos el numero de coordinacion por sitio entre 1 y   \n";
 std::cout << "      100 cada 10 pasos.                                                       \n";
}

void SiteCoord::Evaluate(Configuration & con, Potential & pot)
{
 assert(&pot != 0);
 BasicParticleSet & atoms = con.Atoms();
 const long int natoms = atoms.Size();
 Matrix & m = CurrentValue();
 m = Matrix(3, natoms);
 m.SetLabel(0, "Atom");
 m.SetLabel(1, "Neighbors");
 m.SetLabel(2, "Unique bonds");

 NeighborList * neighbormatrix = new NeighborList[natoms];
 DebugStream() << "-> Computing coordination numbers per site, over " << natoms << " atoms\n";
 DebugStream() << "-> Building neighbor lists\n";
 for (long i=0;i<natoms;++i) neighbormatrix[i] = con.Neighbors(i, true, rcut);

 DebugStream() << "-> Building index table\n";
 // Construct an "index table" so we don't have to depend on Atom::Index()
 std::map<BasicAtom *, long int> indices;
 for (long int i=0;i<atoms.Size();++i) indices[&atoms[i]] = i;
 // 

 for (long int i=0;i<natoms;++i)
 {
  long int nneigh = 0, ubonds = 0;
  NeighborList & nlist = neighbormatrix[i];
  for (long int it=0;it<nlist.Size();++it)
  {
   const AtomPair & nn = nlist[it];
   nneigh++;
   if (indices[nn.j] > i) ubonds++;
  }
  m.Set(0, i, (double)i);
  m.Set(1, i, (double)nneigh);
  m.Set(2, i, (double)ubonds);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SiteCoord(args); }
void destroy(Plugin * m) { delete m; }

