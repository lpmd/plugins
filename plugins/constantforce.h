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
  // Constructor y Destructor
  ConstantForcePotential(std::string args);
  ~ConstantForcePotential();

  double energy(lpmd::SimulationCell & sc);
  void UpdateForces(lpmd::SimulationCell & sc);

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;

 private:
  double fx, fy, fz;

};


#endif

