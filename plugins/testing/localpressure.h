//
//
//

#ifndef __LOCPRESS_H__
#define __LOCPRESS_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class LocalPressure: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  LocalPressure(std::string args);
  ~LocalPressure();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo localpressure
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    int n[3];
    double rcut;
    lpmd::Matrix * m;
    bool do_average;
};

#endif

