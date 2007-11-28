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
   VelocityVerlet(std::string args);
   ~VelocityVerlet();

   // 
   void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);

   void AdvancePositions(lpmd::SimulationCell & sc);
   void AdvanceVelocities(lpmd::SimulationCell & sc);

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;
};

#endif


