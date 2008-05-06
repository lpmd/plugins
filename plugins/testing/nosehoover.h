//
//
//

#ifndef __NOSEHOOVER_H__
#define __NOSEHOOVER_H__

#include <lpmd/vector.h>
#include <lpmd/twostepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

#include <vector>

class NoseHoover: public lpmd::TwoStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
   //Metodos Generales
   NoseHoover(std::string args);
   ~NoseHoover();
   void ShowHelp() const;
   std::string Keywords() const;

   //Metodos Propios del Modulo NoseHoover
   void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);
   void AdvancePosition(lpmd::SimulationCell & sc, long i);
   void AdvanceVelocity(lpmd::SimulationCell & sc, long i);

 private:
      double q;
      double temp;
      double friction;
      std::vector<lpmd::Vector> auxlist;
};

#endif


