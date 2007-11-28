//
//
//

#ifndef __EULER_H__
#define __EULER_H__

#include <lpmd/onestepintegrator.h>
#include <lpmd/applicable.h>
#include <lpmd/plugin.h>

class Euler: public lpmd::OneStepIntegrator, public lpmd::IApplicable, public lpmd::Module
{
 public:
   // 
   Euler(std::string args);
   ~Euler();

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   void Advance(lpmd::SimulationCell & sc);
};

#endif


