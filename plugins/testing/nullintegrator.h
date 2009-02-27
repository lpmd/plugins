//
//
//

#ifndef __NULLINTEGRATOR_H__
#define __NULLINTEGRATOR_H__

#include <lpmd/integrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class NullIntegrator: public lpmd::Integrator, public lpmd::IApplicable, public lpmd::Module
{
 public: 
  //Metodos Generales
  NullIntegrator(std::string args);
  ~NullIntegrator();
  void ShowHelp() const;

  //Metodos propios de modulo nullintegrator
  void Advance(lpmd::SimulationCell & sc, lpmd::Potential & p);

};


#endif


