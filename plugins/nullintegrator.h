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
   NullIntegrator(std::string args);
   ~NullIntegrator();

   //
   void Advance(lpmd::SimulationCell & sc, lpmd::Potential & p);

   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;
};


#endif


