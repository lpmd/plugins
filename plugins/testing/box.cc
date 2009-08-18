//
//
//

#include "box.h"

#include <lpmd/simulation.h>
#include <lpmd/refparticleset.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

class BoxSelector: public Selector<BasicParticleSet>
{
 public:
   BoxSelector(const Box & b): box(b) { }

   const BasicParticleSet & SelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
      if (box.IsInside(ps[i].Position())) innerps.Append(ps[i]);
    return innerps;
   }

   const BasicParticleSet & InverseSelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
      if (!box.IsInside(ps[i].Position())) innerps.Append(ps[i]);
    return innerps;
   }

 private:
   Box box;
   RefParticleSet innerps;
};

BoxFilter::BoxFilter(std::string args): Plugin("box", "1.0"), selector(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("x","0-0");
 DefineKeyword("y","0-0");
 DefineKeyword("z","0-0");
 DefineKeyword("except", "");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 except = params["except"];
 lpmd::Array<std::string> sx = StringSplit(params["x"],'-');
 lpmd::Array<std::string> sy = StringSplit(params["y"],'-');
 lpmd::Array<std::string> sz = StringSplit(params["z"],'-');
 if (sx.Size()!=2 || sx.Size()!=2 || sx.Size()!=2) throw PluginError("box", "Bad, settings in parameters");
 x[0] = atof(sx[0].c_str());x[1] = atof(sx[1].c_str());
 y[0] = atof(sy[0].c_str());y[1] = atof(sy[1].c_str());
 z[0] = atof(sz[0].c_str());z[1] = atof(sz[1].c_str());
}

BoxFilter::~BoxFilter() { delete selector; }

void BoxFilter::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " filter box x=0-10 y=15-20 z=5-10                                              \n";
 std::cout << '\n';
}

Selector<BasicParticleSet> & BoxFilter::CreateSelector()
{
 Box box_region(x[0], x[1], y[0], y[1], z[0], z[1]);
 if (selector != 0) delete selector;
 selector = new BoxSelector(box_region);
 return *selector;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new BoxFilter(args); }
void destroy(Plugin * m) { delete m; }

