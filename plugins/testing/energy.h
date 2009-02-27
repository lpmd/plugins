//
//
//

#ifndef __ENERGY_PROP_H__
#define __ENERGY_PROP_H__

#include <lpmd/vector.h>
#include <lpmd/scalarvalue.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class Energy: public lpmd::ScalarValue, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Energy(std::string args);
  virtual ~Energy();
  void ShowHelp() const;

  //Metodos Propios del modulo energy
  std::string Provides() const;
  double GetProperty(const std::string & name);
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot); 
  const double & Value() const;

 private:
   double ekin, epot, etot, temp;
   lpmd::Vector pv;
};

#endif


