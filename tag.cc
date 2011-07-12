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
   TagSelector(std::string n,std::string v) {name = n;value =v;}

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
   std::string name, value;
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
 DefineKeyword("except", "");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 name = params["name"];
 value = params["value"];
 except = params["except"];
}

TagFilter::~TagFilter() { delete selector; }

void TagFilter::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = tag                                                      \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to identify atoms with an specific Tag.              \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      name          : The name of the tag that you want to compare.            \n";
 std::cout << "      value         : Tag assignation status True/False.                       \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " apply color over tag name=floor value=true                                    \n";
 std::cout << "      In this case we set a 'color' (or anything) to the atoms with the tag    \n";
 std::cout << "      floor in a true status.                                                  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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

