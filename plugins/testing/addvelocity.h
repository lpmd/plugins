//
//
//

#ifndef __EXTRAVEL_H__
#define __EXTRAVEL_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class ExtraVelModifier: public lpmd::SystemModifier, public lpmd::Plugin
{
 public:
  //Metodos Generales
  ExtraVelModifier(std::string args);
  ~ExtraVelModifier();
  void ShowHelp() const;

  void Apply(lpmd::Simulation & sim);

 private:
  lpmd::Vector velocity;
};

#endif



