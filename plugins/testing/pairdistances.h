//
//
//

#ifndef __PD_H__
#define __PD_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class PairDistances: public StoredValue<Matrix>, public InstantProperty, public Module
{
 public:
  //Metodos Generales
  PairDistances(std::string args);
  void ShowHelp() const;

  //Metodos Propios de modulo pairdistances
  void Evaluate(Configuration & conf, Potential & pot);

 private:
    double rcut;
};

#endif

