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

using namespace lpmd;

class Beeman: public lpmd::TwoStepIntegrator, public lpmd::Stepper, public lpmd::Module
{
 public:
   //Metodos Generales
   Beeman(std::string args);
   ~Beeman();
   void ShowHelp() const;

   //Metodos Propios del Modulo Beeman
   void Initialize(Simulation & sim, Potential & p);
   void AdvancePosition(Simulation & sim, long i);
   void AdvanceVelocity(Simulation & sim, long i);

 private:
      std::vector<lpmd::Vector> auxlist;
};

#endif


