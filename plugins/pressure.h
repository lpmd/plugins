//
//
//

#ifndef __PRESSURE_PROP_H__
#define __PRESSURE_PROP_H__

#include <lpmd/scalarvalue.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class Pressure: public lpmd::ScalarValue, public lpmd::InstantProperty, public lpmd::Module
{
 public:
   Pressure(std::string args);
   virtual ~Pressure();

   // From Module
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   //
   void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);
   const double & Value() const;

 private:
   double press;
};

#endif


