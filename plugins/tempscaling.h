//
//
//

#ifndef __TEMPSCALING_H__
#define __TEMPSCALING_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class TempScalingModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
   TempScalingModifier(std::string args);
   ~TempScalingModifier();

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   void Apply(lpmd::MD & md);

  private:
    double fromtemp;
    double totemp;
};

#endif



