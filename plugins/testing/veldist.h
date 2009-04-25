//
//
//

#ifndef __VELDIST_H__
#define __VELDIST_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class VelDist: public lpmd::Value<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  VelDist(std::string args);
  ~VelDist();
  void ShowHelp() const;

  //Metodos Propios de modulo veldist
  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int bins;
};

#endif

