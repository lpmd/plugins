//
//
//

#ifndef __VERLET_H__
#define __VERLET_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/stepper.h>
#include <lpmd/plugin.h>

class Verlet: public lpmd::OneStepIntegrator, public lpmd::Stepper, public lpmd::Module
{
 public:
  //Metodos Generales
  Verlet(std::string args);
  ~Verlet();
  void ShowHelp() const;

  //Metodos propios modulo verlet
  void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);
  void Advance(lpmd::SimulationCell & sc, long i);
};


#endif


