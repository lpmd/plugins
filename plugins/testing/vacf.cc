//
//
//

#include "vacf.h"
#include <lpmd/properties.h>

using namespace lpmd;

Vacf::Vacf(std::string args): Module("vacf")
{
 ParamList & param = (*this);
 m = NULL;
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 AssignParameter("dt","0");
 //
 ProcessArguments(args);
 dt = param["dt"];
}

Vacf::~Vacf() { if (m != NULL) delete m; }

void Vacf::Evaluate(lpmd::SimulationHistory & hist, Potential & pot)
{
 if(m!=NULL) delete []m;
 m=vacf(hist,pot,dt);
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Vacf(args); }
void destroy(Module * m) { delete m; }
