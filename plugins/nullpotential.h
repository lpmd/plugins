//
//
//

#ifndef __NULLPOTENTIAL_H__
#define __NULLPOTENTIAL_H__

#include <lpmd/potential.h>
#include <lpmd/plugin.h>

class NullPotential: public lpmd::Potential, public lpmd::Module
{
 public:
   // Constructor y Destructor
   NullPotential(std::string args); 
   ~NullPotential();

  double energy(lpmd::SimulationCell & sc);
  void UpdateForces(lpmd::SimulationCell & sc);

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;

};


#endif

