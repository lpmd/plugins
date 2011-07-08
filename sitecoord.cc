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
 // hasta aqui los valores por omision
 ProcessArguments(args);
 rcut = params["rcut"];
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
}

SiteCoord::~SiteCoord() { }

void SiteCoord::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = sitecoord                                                \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to evaluate the coordination number by site in a     \n";
 std::cout << "      simulation cell.                                                         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      rcut          : Set the maximum cutoff for the evaluation.               \n";
 std::cout << "      output        : output filename.                                         \n";
 std::cout << "      average       : True or False in case that realize an average over the   \n";
 std::cout << "                      configurations.                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " use sitecoord                                                                 \n";
 std::cout << "     output sitecoord.dat                                                      \n";
 std::cout << "     rcut 4.0                                                                  \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #Applying the plugin :                                                        \n";  
 std::cout << " property sitecoord start=1 each=10 end=100                                    \n";
 std::cout << "      With this we calculate the coordination number by site between 1 and 100 \n";
 std::cout << "      each 10 timestep.                                                        \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 delete [] neighbormatrix;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SiteCoord(args); }
void destroy(Plugin * m) { delete m; }

