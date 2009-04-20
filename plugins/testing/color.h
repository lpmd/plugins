//
//
//

#ifndef __COLOR_H__
#define __COLOR_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class ColorModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  ColorModifier(std::string args);
  ~ColorModifier();
  void ShowHelp() const;

  //Metodos Propios
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

 private:
  double fromtemp, totemp;
};

#endif



