//
//
//

#ifndef __PD_H__
#define __PD_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class PairDistances: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  PairDistances(std::string args);
  ~PairDistances();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo pairdistances
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    double rcut;
};

#endif

