//
//
//

#ifndef __EXTRAVEL_H__
#define __EXTRAVEL_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class ExtraVelModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  ExtraVelModifier(std::string args);
  ~ExtraVelModifier();
  void ShowHelp() const;

  //Metodos Propios
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);
};

#endif



