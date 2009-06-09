//
//
//

#ifndef __METROPOLI_H__
#define __METROPOLI_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/stepper.h>
#include <lpmd/plugin.h>

class Metropoli: public lpmd::OneStepIntegrator, public lpmd::Plugin
{
 public:
  //Metodos Generales
  Metropoli(std::string args);
  ~Metropoli();
  void ShowHelp() const;

  //Metodos propios modulo verlet
  void Initialize(lpmd::Simulation & sim, lpmd::Potential & p);
  void Advance(lpmd::Simulation & sim, long i);
 private:
  double Temp;
  double rand;
};


#endif


