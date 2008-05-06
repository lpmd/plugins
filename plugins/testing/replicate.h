//
//
//

#ifndef __REPLICATE_H__
#define __REPLICATE_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class ReplicateModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:

  //Metodos Generales 
  ReplicateModifier(std::string args);
  ~ReplicateModifier();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo cellscaling
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

 private:
  unsigned long nx, ny, nz;
};

#endif

