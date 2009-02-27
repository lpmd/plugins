//
//
//

#ifndef __SELECT_ATOMS_H__
#define __SELECT_ATOMS_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>
#include <lpmd/vector.h>

class SelectAtomsModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:

  //Metodos Generales 
  SelectAtomsModifier(std::string args);
  ~SelectAtomsModifier();
  void SetParameter(std::string name);
  void Show(std::ostream & os) const;
  void ShowHelp() const;

  //Metodos Propios de modulo selectatoms
  void Apply(lpmd::SimulationCell & sc);
  void Apply(lpmd::MD & md);

  private:
    std::string mode;
    lpmd::Vector vmin, vmax, center, vel, acc;
    double radius;
    long p0, p1;
    bool outside;
};

#endif

