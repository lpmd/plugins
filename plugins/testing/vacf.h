//
//
//

#ifndef __VACF_H__
#define __VACF_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/temporalproperty.h>
#include <lpmd/plugin.h>

class Vacf: public lpmd::Value<lpmd::Matrix>, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Vacf(std::string args);
  ~Vacf();

  //Metodos propios de modulo vacf
  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(const std::vector<lpmd::SimulationCell> & simcell, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
  double dt;
};

#endif

