//
//
//

#ifndef __CELLSCALING_H__
#define __CELLSCALING_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class CellScalingModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:

  //Metodos Generales 
  CellScalingModifier(std::string args);
  ~CellScalingModifier();
  void Show(std::ostream & os) const;
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo cellscaling
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

 private:
  int axis;
  double percent;
};

#endif



