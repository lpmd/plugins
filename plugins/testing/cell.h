//
//
//

#ifndef __CELL_PROP_H__
#define __CELL_PROP_H__

#include <lpmd/scalarvalue.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class CellProp: public lpmd::ScalarValue, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  CellProp(std::string args);
  virtual ~CellProp();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios del modulo cell
  std::string Provides() const;
  double GetProperty(const std::string & name);
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot); 
  const double & Value() const;

 private:
   double volume, a, b, c, dens, partdens;
   long nat;
};

#endif


