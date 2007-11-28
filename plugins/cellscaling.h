//
//
//

#ifndef __CELLSCALING_H__
#define __CELLSCALING_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/plugin.h>

class CellScalingModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:

   CellScalingModifier(std::string args);
   ~CellScalingModifier();

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   void Apply(lpmd::MD & md);

  private:
    int axis;
    double percent;
};

#endif



