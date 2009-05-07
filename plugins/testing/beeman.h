//
//
//

#ifndef __BEEMAN_H__
#define __BEEMAN_H__

#include <lpmd/vector.h>
#include <lpmd/twostepintegrator.h>
#include <lpmd/stepper.h>
#include <lpmd/plugin.h>

#include <vector>

class Beeman: public lpmd::TwoStepIntegrator, public lpmd::Stepper, public lpmd::Module
{
 public:
   //Metodos Generales
   Beeman(std::string args);
   ~Beeman();
   void ShowHelp() const;

   //Metodos Propios del Modulo Beeman
   void Initialize(lpmd::SimulationCell & sc, lpmd::Potential & p);
   void AdvancePosition(lpmd::SimulationCell & sc, long i);
   void AdvanceVelocity(lpmd::SimulationCell & sc, long i);

 private:
      std::vector<lpmd::Vector> auxlist;
};

#endif


