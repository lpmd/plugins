//
//
//

#ifndef __CORDNUM_H__
#define __CORDNUM_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

#include <map>

class CordNum: public lpmd::Value<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  CordNum(std::string args);
  ~CordNum();
  void SetParameter(std::string name);
  void Show(std::ostream & os) const;
  void ShowHelp() const;

  //Metodos Propios del Modulo cordnum
  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
  int nb;
  std::map<std::string, double> rcut;
  int na;
  std::vector<std::string> satoms; //atomic symbols.
  bool do_average;
};

#endif

