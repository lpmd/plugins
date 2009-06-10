//
//
//

#ifndef __UNDOPBC_H__
#define __UNDOPBC_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>
#include <lpmd/vector.h>

class UndoPBCModifier: public lpmd::SystemModifier, public lpmd::Plugin
{
 public:

  //Metodos Generales 
  UndoPBCModifier(std::string args);
  ~UndoPBCModifier();
  void Show(std::ostream & os) const;
  void ShowHelp() const;

  //Metodos Propios de modulo cellscaling
  void Apply(lpmd::Simulation & sim);

 private:
  int counter;
  lpmd::Configuration old;
};

#endif



