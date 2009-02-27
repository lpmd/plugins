//
//
//

#ifndef __GDR_H__
#define __GDR_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class Gdr: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Gdr(std::string args);
  ~Gdr();
  void ShowHelp() const;

  //Metodos Propios de modulo gdr
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    double rcut;
    int nb;
    bool do_average;
};

#endif

