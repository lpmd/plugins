//
//
//

#ifndef __TEMPSCALING_H__
#define __TEMPSCALING_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>
#include <lpmd/vector.h>

class ThermalNeedleModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  ThermalNeedleModifier(std::string args);
  ~ThermalNeedleModifier();
  void ShowHelp() const;

  //Metodos Propios
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

 private:
  double temperature;
  lpmd::Vector center;
  double radius;
};

#endif



