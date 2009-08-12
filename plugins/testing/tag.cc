//
//
//

#include "tag.h"

#include <lpmd/simulation.h>
#include <lpmd/refparticleset.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

class TagSelector: public Selector<BasicParticleSet>
{
 public:
   TagSelector(std::string n,bool v) {name = n;value =v;}

   const BasicParticleSet & SelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i)
    { 
     if (ps.Have(ps[i],Tag(name)) && (ps.GetTag(ps[i], Tag(name))==value)) innerps.Append(ps[i]);
    }
    return innerps;
   }

   const BasicParticleSet & InverseSelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
    {
     if (!ps.Have(ps[i],Tag(name)) || !(ps.GetTag(ps[i], Tag(name))==value)) innerps.Append(ps[i]);
    }
    return innerps;
   }

 private:
   std::string name;
   bool value;
   RefParticleSet innerps;
};

TagFilter::TagFilter(std::string args): Plugin("tag", "1.0"), selector(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("name","fixedvel");
 DefineKeyword("value","true");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 name = params["name"];
 value = bool(params["value"]);
}

TagFilter::~TagFilter() { delete selector; }

void TagFilter::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " filter tag name=fixedpos value=true                                          \n";
 std::cout << '\n';
}

Selector<BasicParticleSet> & TagFilter::CreateSelector()
{
 if (selector != 0) delete selector;
 selector = new TagSelector(name,value);
 return *selector;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new TagFilter(args); }
void destroy(Plugin * m) { delete m; }

