//
//
//

#ifndef __LEAPFROG_H__
#define __LEAPFROG_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class Leapfrog: public lpmd::OneStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
  //Metodos Generales
  Leapfrog(std::string args);
  ~Leapfrog();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos propios modulo verlet
  void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);
  void Advance(lpmd::SimulationCell & sc, long i);
};


#endif


