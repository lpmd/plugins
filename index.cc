//
//
//

#include "index.h"

#include <lpmd/simulation.h>
#include <lpmd/refparticleset.h>
#include <lpmd/util.h>

#include <iostream>
#include <string>

using namespace lpmd;

class IndexSelector: public Selector<BasicParticleSet>
{
 public:
   IndexSelector(lpmd::Array<int> i) {index = i;}

   const BasicParticleSet & SelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i)
    {
     for (long int j=0;j<index.Size();++j)
     {
      if (i==index[j]) innerps.Append(ps[i]);
     }
    } 
    return innerps;
   }

   const BasicParticleSet & InverseSelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i)
    {
     for (long int j=0;j<index.Size();++j)
     {
      if (i!=index[j]) innerps.Append(ps[i]);
     }
    } 
    return innerps;
   }

 private:
   lpmd::Array<int> index;
   RefParticleSet innerps;
};

IndexFilter::IndexFilter(std::string args): Plugin("index", "1.0"), selector(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("index","0-1");
 DefineKeyword("except", "");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 std::string tmp = params["index"];
 except = params["except"];
 size_t found = std::string::npos - 1;
 int ct=-1;
 while(found!=std::string::npos)
 {
  ct++;
  found=tmp.find(",");
 }
 if (ct > 0) index = StringSplit(params["index"],',');
 else if (ct == 0)
 {
  Array<std::string> limits = StringSplit(params["index"],'-');
  if (limits.Size() != 2) throw PluginError("index", "Wrong specification of \"index\" range");
  for (int i=atoi(limits[0].c_str());i<=atoi(limits[1].c_str());++i)
      index.Append(ToString<int>(i));
 }
 if (args != "dummyargument") 
    DebugStream() << "-> Reading " << index.Size() << " indices to filter." << '\n';
}

IndexFilter::~IndexFilter() { delete selector; }

void IndexFilter::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = index                                                    \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This module is used to select atoms by their index in the list of atoms  \n";
 std::cout << "      of the cell. It can be called with the 'filter' or 'over' keyword.       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      index         : Sets the indices (separated by comma) or range (separated\n";
 std::cout << "                      by a minus (-) sign) in the atoms list that are wanted   \n";
 std::cout << "                      to be considered. The index of an atom corresponds to the\n";
 std::cout << "                      place in which the atom appears in the input file (the   \n";
 std::cout << "                      first lines will corespond to the atoms 1,2,3,4,etc.).   \n";
 std::cout << "      start         : Determines in which step the plugin begins to be applied.\n";
 std::cout << "      end           : Determines in which step the plugin ceases to be applied.\n";
 std::cout << "      each          : Determines how often (each how many time-steps) the      \n";
 std::cout << "                      plugin must be applied.                                  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " #Applying plugin :                                                            \n";
 std::cout << " filter index index=5,6,8,10                                                   \n";
 std::cout << " apply myplugin over index index=5-10 start=0 end=1 each=1                   \n\n";
 std::cout << "      The plugin is used to eliminate (filter) four atoms: the 5th, 6th, 8th   \n";
 std::cout << "      and 10th of the atoms list in the first case, and to apply the 'myplugin'\n";
 std::cout << "      plugin to all atoms with an index in the range from 5 to 10 in the second\n";
 std::cout << "      case ('myplugin' must be loaded previously).                             \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

Selector<BasicParticleSet> & IndexFilter::CreateSelector()
{
 if (selector != 0) delete selector;
 lpmd::Array<int> in;
 in.Clear();
 for(int i=0;i<index.Size();++i)
 {
  in.AppendUnique(atoi(index[i].c_str()));
 }
 selector = new IndexSelector(in);
 return *selector;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new IndexFilter(args); }
void destroy(Plugin * m) { delete m; }

