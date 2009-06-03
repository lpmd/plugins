//
//
//

#ifndef __CORDNUMFUNC_H__
#define __CORDNUMFUNC_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>

class CordNumFunc: public lpmd::StoredValue<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  CordNumFunc(std::string args);
  ~CordNumFunc();
  void SetParameter(std::string name);
  void Show(std::ostream & os) const;
  void ShowHelp() const;

  //Metodos Propios de la clase cordnumfunc
  void Evaluate(lpmd::Configuration & con, lpmd::Potential & pot);

 private:
  int nb;
  double cut;
  int na;
  std::vector<std::string> satoms; //atomic symbols.
  bool do_average;
};

#endif

