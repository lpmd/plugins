//
//
//

#ifndef __DISPLACE_SM_H__
#define __DISPLACE_SM_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/vector.h>
#include <lpmd/plugin.h>

class DisplaceModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  DisplaceModifier(std::string args);
  ~DisplaceModifier();

  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);
 
 private:
  lpmd::Vector offset;
};

#endif

