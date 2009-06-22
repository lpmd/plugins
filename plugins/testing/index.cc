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

IndexFilter::IndexFilter(std::string args): Plugin("element", "1.0"), selector(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("index","0-1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 std::string tmp = params["index"];
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
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " filter index index=5,6,8,10                                                   \n";
 std::cout << " filter index index=5-10                                                       \n";
 std::cout << '\n';
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

