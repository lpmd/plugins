//
//
//

#ifndef __PD_H__
#define __PD_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class PairDistances: public Value<Matrix>, public InstantProperty, public Module
{
 public:
  //Metodos Generales
  PairDistances(std::string args);
  ~PairDistances();
  void ShowHelp() const;

  //Metodos Propios de modulo pairdistances
  const Matrix & CurrentValue() const { return *m; }
  void Evaluate(Configuration & conf, Potential & pot);

 private:
    Matrix * m;
    double rcut;
};

#endif

