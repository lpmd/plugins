//
//
//

#ifndef __EULER_H__
#define __EULER_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class Euler: public lpmd::OneStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
  //Metodos Generales
  Euler(std::string args);
  ~Euler();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios del modulo euler
  void Advance(lpmd::SimulationCell & sc, long i);
};

#endif


