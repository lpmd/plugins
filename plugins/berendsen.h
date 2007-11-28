//
//
//

#ifndef __BERENDSEN_H__
#define __BERENDSEN_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class BerendsenModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:

   BerendsenModifier(std::string args);
   ~BerendsenModifier();

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   void Apply(lpmd::MD & md);

  private:
    long stop_thermostat, old_step;
    double fromtemp, totemp, tau;
};

#endif



