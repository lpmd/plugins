//
//
//

#ifndef __VERLET_H__
#define __VERLET_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/stepper.h>
#include <lpmd/plugin.h>

class HardSpheres: public lpmd::OneStepIntegrator, public lpmd::Plugin
{
 public:
  //Metodos Generales
  HardSpheres(std::string args);
  ~HardSpheres();
  void ShowHelp() const;

  //Metodos propios modulo verlet
  void Initialize(lpmd::Simulation & sim, lpmd::Potential & p);
  void Advance(lpmd::Simulation & sim, long i);
};


#endif


