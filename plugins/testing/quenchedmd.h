//
//
//

#ifndef __QUENCHEDMD_H__
#define __QUENCHEDMD_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class QuenchedMDModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  QuenchedMDModifier(std::string args);
  ~QuenchedMDModifier();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios
  void Apply(lpmd::Simulation & sim);
};

#endif



