//
//
//

#ifndef __VERLET_H__
#define __VERLET_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class Verlet: public lpmd::OneStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
   Verlet(std::string args);
   ~Verlet();

   // Overloaded from Integrator
   void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);

   // Overloaded from OneStepIntegrator
   void Advance(lpmd::SimulationCell & sc);

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;
};


#endif


