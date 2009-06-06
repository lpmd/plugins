//
//
//

#include "monitor.h"

#include <lpmd/simulation.h>

#include <iostream>
#include <fstream>

using namespace lpmd;

MonitorVisualizer::MonitorVisualizer(std::string args): Plugin("monitor", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("properties", "step");
 DefineKeyword("output", "-");
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 if (params["output"] != "-") output_stream = new std::ofstream(params["output"].c_str());
 else output_stream = &(std::cout);
 property_array = StringSplit(params["properties"], ',');
}

MonitorVisualizer::~MonitorVisualizer()
{ 
 ParamList & params = (*this);
 if (params["output"] != "-") delete output_stream;
}

void MonitorVisualizer::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " monitor step,temperature start=1 end=1000 each=50                           \n\n";
}

void MonitorVisualizer::Apply(const Simulation & sim)
{ 
 for (int p=0;p<property_array.Size();++p) 
     (*output_stream) << sim.GetTag(sim, Tag(property_array[p])) << "  ";
 (*output_stream) << '\n';
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new MonitorVisualizer(args); }
void destroy(Plugin * m) { delete m; }

