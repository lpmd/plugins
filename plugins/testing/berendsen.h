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
   //Metodos Generales
   BerendsenModifier(std::string args);
   ~BerendsenModifier();
   void ShowHelp() const;

   //Metodos Propios del Modulo Berendsen
   void Apply(lpmd::SimulationCell & sc);
   void Apply(lpmd::MD & md);

  private:
    long stop_thermostat, old_step;
    double fromtemp, totemp, tau;
};

#endif



