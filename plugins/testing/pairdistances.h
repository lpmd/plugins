//
//
//

#ifndef __PD_H__
#define __PD_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class PairDistances: public lpmd::Value<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  PairDistances(std::string args);
  ~PairDistances();
  void ShowHelp() const;

  //Metodos Propios de modulo pairdistances
  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    double rcut;
};

#endif

