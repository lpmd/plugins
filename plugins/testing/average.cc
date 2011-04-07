//
//
//

#include "average.h"

#include <lpmd/simulation.h>

#include <iostream>
#include <fstream>

using namespace lpmd;

AverageVisualizer::AverageVisualizer(std::string args): Plugin("average", "1.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("properties", "step");
 DefineKeyword("output", "-");
 DefineKeyword("interval", "100");
 ProcessArguments(args);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 if (params["output"] != "-") output_stream = new std::ofstream(params["output"].c_str());
 else output_stream = &(std::cout);
 interval = int(params["interval"]);
 property_array = StringSplit(params["properties"], ',');
 property_averages = new double[property_array.Size()];
 last_value = new double[property_array.Size()];
 for (int i=0;i<property_array.Size();++i) property_averages[i] = last_value[i] = 0.0;
 average_count = 0;
}

AverageVisualizer::~AverageVisualizer()
{ 
 ParamList & params = (*this);
 if (params["output"] != "-") delete output_stream;
 delete [] property_averages;
 delete [] last_value;
}

void AverageVisualizer::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      The plugin monitor the average of specific properties during the         \n";
 std::cout << " simulation process.                                                         \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "                     Availables properties :                                   \n";
 std::cout << "                                           step,temperature,volume,pressure,   \n";
 std::cout << "                                           volume-per-atom,cell-a,cell-b,cell-c\n";
 std::cout << "                                           particle-density,density,momentum   \n";
 std::cout << "                                           px,py,pz,potential-energy,          \n";
 std::cout << "                                           kinetic-energy,total-energy,        \n";
 std::cout << "                                           virial-pressure,kinetic-pressure    \n";
 std::cout << "                                           sxx,sxy,sxz,etc.                    \n";
 std::cout << "                     For a detailed information about units, refeer to the     \n";
 std::cout << "                     manual of lpmd in General Proerties section.              \n";
 std::cout << '\n';
 std::cout << " Example           >>                                                          \n";
 std::cout << " average temperature,sxx interval=100 start=0 end=-1 each=50 output=out.dat  \n\n";
}

void AverageVisualizer::Apply(const Simulation & sim)
{ 
 for (int p=0;p<property_array.Size();++p) 
 {
  double value = double(Parameter(sim.GetTag(sim, Tag(property_array[p]))));
  if (average_count < interval)
  {
   property_averages[p] += value;
  }
  else
  {
   property_averages[p] -= last_value[p];
   property_averages[p] += value;
  }
  (*output_stream) << property_averages[p]/double(average_count) << " ";
  last_value[p] = value;
 }
 (*output_stream) << '\n';
 if (average_count < interval) average_count++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new AverageVisualizer(args); }
void destroy(Plugin * m) { delete m; }

