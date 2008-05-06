//
//
//

#ifndef __CORDNUMFUNC_H__
#define __CORDNUMFUNC_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class CordNumFunc: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  CordNumFunc(std::string args);
  ~CordNumFunc();
  void SetParameter(std::string name);
  void Show(std::ostream & os) const;
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de la clase cordnumfunc
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
  int nb;
  double cut;
  int na;
  std::vector<std::string> satoms; //atomic symbols.
  bool do_average;
};

#endif

