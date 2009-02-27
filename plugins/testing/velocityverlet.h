//
//
//

#ifndef __VELOCITYVERLET_H__
#define __VELOCITYVERLET_H__

#include <lpmd/twostepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class VelocityVerlet: public lpmd::TwoStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
  //Metodos Generales
  VelocityVerlet(std::string args);
  ~VelocityVerlet();
  void ShowHelp() const;

  //Metodos Propios modulo velocityverlet
  void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);
  void AdvancePosition(lpmd::SimulationCell & sc, long i);
  void AdvanceVelocity(lpmd::SimulationCell & sc, long i);
};

#endif
