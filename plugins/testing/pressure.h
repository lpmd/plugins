//
//
//

#ifndef __PRESSURE_PROP_H__
#define __PRESSURE_PROP_H__

#include <lpmd/value.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class Pressure: public lpmd::Value<double>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Pressure(std::string args);
  virtual ~Pressure();
  void ShowHelp() const;

  //Metodos Propios del modulo pressure
  std::string Provides() const;
  double GetProperty(const std::string & name);
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot); 
  const double & CurrentValue() const;

 private:
   double press, kpress, vpress;
   double s[3][3];
};

#endif


