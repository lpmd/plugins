//
//
//

#ifndef __DENSITYPROFILE_H__
#define __DENSITYPROFILE_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class DensityProfile: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  DensityProfile(std::string args);
  ~DensityProfile();
  void SetParameter(std::string name);
  void Show(std::ostream & os) const;
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo gdr
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int bins;
    int axis;
    double range[3][2];
    bool do_average;
};

#endif

