//
//
//

#include "nullmixer.h"

#include <iostream>

using namespace lpmd;

NullMixer::NullMixer(std::string args): Module("nullmixer")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 ProcessArguments(args);
}

NullMixer::~NullMixer() { }

void NullMixer::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << "  Aplicando el Modulo                                                          \n";
 std::cout << " apply nullmixer                                                               \n";
}

Configuration & NullMixer::Apply(Configuration & config1, Configuration & config2)
{
 return config1;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new NullMixer(args); }
void destroy(Module * m) { delete m; }

