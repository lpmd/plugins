//
//
//

#include "sphere.h"

#include <lpmd/simulation.h>
#include <lpmd/refparticleset.h>
#include <lpmd/util.h>

#include <iostream>

using namespace lpmd;

class SphereSelector: public Selector<BasicParticleSet>
{
 public:
   SphereSelector(const Sphere & sph): sphere(sph) { }

   const BasicParticleSet & SelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
      if (sphere.IsInside(ps[i].Position())) innerps.Append(ps[i]);
    return innerps;
   }

   const BasicParticleSet & InverseSelectFrom(const BasicParticleSet & ps) 
   { 
    innerps.Clear();
    for (long int i=0;i<ps.Size();++i) 
      if (!sphere.IsInside(ps[i].Position())) innerps.Append(ps[i]);
    return innerps;
   }

 private:
   Sphere sphere;
   RefParticleSet innerps;
};

SphereFilter::SphereFilter(std::string args): Plugin("sphere", "1.0"), selector(0)
{
 ParamList & params = (*this);
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("radius");
 DefineKeyword("center");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 radius = double(params["radius"]);
 center = Vector(params["center"].c_str());
}

SphereFilter::~SphereFilter() { delete selector; }

void SphereFilter::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " filter sphere radius=5.0 center=10.0,10.0,10.0                                \n";
 std::cout << '\n';
}

Selector<BasicParticleSet> & SphereFilter::CreateSelector()
{
 Sphere spheric_region(center, radius);
 if (selector != 0) delete selector;
 selector = new SphereSelector(spheric_region);
 return *selector;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new SphereFilter(args); }
void destroy(Plugin * m) { delete m; }

