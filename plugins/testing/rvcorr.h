//
//
//

#ifndef __RVCORR_H__
#define __RVCORR_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class RVCorr: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  RVCorr(std::string args);
  ~RVCorr();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo rvcorr
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    double rcut;
    int nb;
    bool do_average;
};

#endif

