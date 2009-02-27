//
//
//

#ifndef __CONSTANT_FORCE_POTENTIAL_H__
#define __CONSTANT_FORCE_POTENTIAL_H__

#include <lpmd/potential.h>
#include <lpmd/plugin.h>

class ConstantForcePotential: public lpmd::Potential, public lpmd::Module
{
 public:
  //Metodos Generales
  ConstantForcePotential(std::string args);
  ~ConstantForcePotential();
  void SetParameter(std::string name);
  void ShowHelp() const;

  //Metodos Propios del modulo constantforce
  double energy(lpmd::SimulationCell & sc);
  void UpdateForces(lpmd::SimulationCell & sc);

 private:
  double fx, fy, fz;

};


#endif

