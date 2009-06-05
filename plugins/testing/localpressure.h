//
//
//

#ifndef __LOCPRESS_H__
#define __LOCPRESS_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>
#include <lpmd/simulation.h>

class LocalPressure: public lpmd::StoredValue<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  LocalPressure(std::string args);
  ~LocalPressure();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo localpressure
  void Evaluate(lpmd::Configuration & con, lpmd::Potential & pot);

 private:
    int n[3];
    double rcut;
    bool do_average;
};

#endif

