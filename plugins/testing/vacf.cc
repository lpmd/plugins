//
//
//

#include "vacf.h"
#include "plugincommon.h"

using namespace lpmd;

Vacf::Vacf(std::string args): Module("vacf")
{
 m = NULL;
 AssignParameter("dt","0");
 ProcessArguments(args);
 dt = GetDouble("dt");
}

Vacf::~Vacf() { if (m != NULL) delete m; }

std::string Vacf::Keywords() const { return "dt"; }

void Vacf::Evaluate(const std::vector<SimulationCell> & simcell, Potential & pot)
{
 if(m!=NULL) delete[]m;
 m=vacf(simcell,pot,dt);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Vacf(args); }
void destroy(Module * m) { delete m; }
