//
//
//

#include "box.h"

#include <lpmd/simulation.h>
#include <lpmd/refparticleset.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

#define DELTA (1.0E+12)

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
 if (params["x"] == "all")
 {
  x[0] = -DELTA;
  x[1] = +DELTA;
 }
 else
 {
  lpmd::Array<std::string> sx = StringSplit(params["x"],'-');
  if (sx.Size() != 2) throw PluginError("box", "Invalid x range");
  x[0] = atof(sx[0].c_str());
  x[1] = atof(sx[1].c_str());
 }

 if (params["y"] == "all")
 {
  y[0] = -DELTA;
  y[1] = +DELTA;
 }
 else
 {
  lpmd::Array<std::string> sy = StringSplit(params["y"],'-');
  if (sy.Size() != 2) throw PluginError("box", "Invalid y range");
  y[0] = atof(sy[0].c_str());
  y[1] = atof(sy[1].c_str());
 }

 if (params["z"] == "all")
 {
  z[0] = -DELTA;
  z[1] = +DELTA;
 } 
 else
 {
  lpmd::Array<std::string> sz = StringSplit(params["z"],'-');
  if (sz.Size() != 2) throw PluginError("box", "Invalid z range");
  z[0] = atof(sz[0].c_str());
  z[1] = atof(sz[1].c_str());
 }
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
 ParamList & params = (*this);
 DebugStream() << "-> Selecting range: ";
 DebugStream() << "X=";
 if (params["x"] != "all") DebugStream() << "[ " << x[0] << ", " << x[1] << " ] ";
 else DebugStream() << "all "; 
 DebugStream() << "Y=";
 if (params["y"] != "all") DebugStream() << "[ " << y[0] << ", " << y[1] << " ] ";
 else DebugStream() << "all "; 
 DebugStream() << "Z=";
 if (params["z"] != "all") DebugStream() << "[ " << z[0] << ", " << z[1] << " ] ";
 else DebugStream() << "all "; 
 DebugStream() << "\n";
 if (selector != 0) delete selector;
 selector = new BoxSelector(box_region);
 return *selector;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new BoxFilter(args); }
void destroy(Plugin * m) { delete m; }

