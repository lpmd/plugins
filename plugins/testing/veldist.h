//
//
//

#ifndef __VELDIST_H__
#define __VELDIST_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class VelDist: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  VelDist(std::string args);
  ~VelDist();
  void ShowHelp() const;

  //Metodos Propios de modulo veldist
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int bins;
};

#endif

