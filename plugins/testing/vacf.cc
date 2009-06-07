//
//
//

#include "vacf.h"
#include <lpmd/properties.h>

using namespace lpmd;

Vacf::Vacf(std::string args): Plugin("vacf", "2.0")
{
 ParamList & param = (*this);
 AssignParameter("dt","0");
 //
 ProcessArguments(args);
 dt = param["dt"];
}

Vacf::~Vacf() { }

void Vacf::Evaluate(lpmd::ConfigurationSet & hist, Potential & pot)
{
 vacf(hist,pot,dt,CurrentValue());
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Vacf(args); }
void destroy(Plugin * m) { delete m; }
