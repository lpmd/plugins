//
//
//

#ifndef __BEEMAN_H__
#define __BEEMAN_H__

#include <lpmd/twostepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class Beeman: public lpmd::TwoStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
   Beeman(std::string args);
   ~Beeman();

   //
   void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);

   void AdvancePositions(lpmd::SimulationCell & sc);
   void AdvanceVelocities(lpmd::SimulationCell & sc);

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   private:
      std::vector<lpmd::Vector> auxlist;
};

#endif


