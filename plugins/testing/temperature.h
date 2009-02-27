//
//
//

#ifndef __TEMPERATURE_SM_H__
#define __TEMPERATURE_SM_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class TemperatureModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  TemperatureModifier(std::string args);
  ~TemperatureModifier();

  void ShowHelp() const;

  //Metodos Propios
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

 private:
  double temp;
};

#endif

