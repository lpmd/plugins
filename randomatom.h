//
//
//

#ifndef __RANDOMATOM_SM_H__
#define __RANDOMATOM_SM_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/vector.h>
#include <lpmd/plugin.h>
#include <lpmd/configuration.h>

using namespace lpmd;

class RandomAtomModifier: public lpmd::SystemModifier, public lpmd::Plugin
{
 public:
  //Metodos Generales
  RandomAtomModifier(std::string args);
  ~RandomAtomModifier();

  void ShowHelp() const;

  //Metodos Propios
  void Apply(Simulation & conf);
 
 private:
  std::string type;
  double value;
  std::string symbol;
  std::string density;
};

#endif

