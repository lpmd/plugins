//
//
//

#ifndef __TEMPSCALING_H__
#define __TEMPSCALING_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class TempScalingModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  TempScalingModifier(std::string args);
  ~TempScalingModifier();
  void ShowHelp() const;

  //Metodos Propios
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

 private:
  double fromtemp, totemp;
};

#endif



