//
//
//

#ifndef __TEMPPROFILE_H__
#define __TEMPPROFILE_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class TempProfile: public lpmd::Value<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  TempProfile(std::string args);
  ~TempProfile();
  void SetParameter(std::string name);
  void Show(std::ostream & os) const;
  void ShowHelp() const;

  //Metodos Propios de modulo gdr
  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int bins;
    int axis;
    double range[3][2];
    bool do_average;
    long int counter;
};

#endif

