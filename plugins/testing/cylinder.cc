//
//
//

#include "cylinder.h"

#include <lpmd/simulation.h>
#include <lpmd/refparticleset.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

class CylinderSelector: public Selector<BasicParticleSet>
{
 public:
   CylinderSelector(const Cylinder & b): cyl(b) { }

   const BasicParticleSet & SelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
      if (cyl.IsInside(ps[i].Position())) innerps.Append(ps[i]);
    return innerps;
   }

   const BasicParticleSet & InverseSelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
      if (!cyl.IsInside(ps[i].Position())) innerps.Append(ps[i]);
    return innerps;
   }

 private:
   Cylinder cyl;
   RefParticleSet innerps;
};

CylinderFilter::CylinderFilter(std::string args): Plugin("cylinder", "1.0"), selector(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("rmax","0");
 DefineKeyword("rmin","0");
 DefineKeyword("origin","<0,0,0>");
 DefineKeyword("endpoint","null");
 DefineKeyword("except", "");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 except = params["except"];
 //if (params["endpoint"] == "null") throw PluginError("cylinder", "Invalid end point");
 //if (params["rmax"] == "0") throw PluginError("rmax", "Invalid maxim radius");
 origin = Vector(params["origin"].c_str());
 S = Vector(params["endpoint"].c_str()) - Vector(params["origin"].c_str());
 rmax = double(params["rmax"]);
 rmin = double(params["rmin"]);
}

CylinderFilter::~CylinderFilter() { delete selector; }

void CylinderFilter::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " filter cylinder origin=<5,5,5> endpoint=<6,6,6> rmax=10 rmin=0                     \n";
 std::cout << '\n';
}

Selector<BasicParticleSet> & CylinderFilter::CreateSelector()
{
 Cylinder cyl_region(S, origin, rmax, rmin);
 ParamList & params = (*this);
 DebugStream() << "-> Selecting cylinder: "<< '\n';
 DebugStream() << "origin   = " << params["origin"] << '\n';
 DebugStream() << "endpoint = " << params["endpoint"] << '\n';
 DebugStream() << "rmax     = " << params["rmax"] << '\n';
 DebugStream() << "rmin     = " << params["rmin"] << '\n';
 DebugStream() << "Volume   = " << cyl_region.Volume() << '\n';
 if (selector != 0) delete selector;
 selector = new CylinderSelector(cyl_region);
 return *selector;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new CylinderFilter(args); }
void destroy(Plugin * m) { delete m; }

