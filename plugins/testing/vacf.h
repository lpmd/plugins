//
//
//

#ifndef __VACF_H__
#define __VACF_H__

#include <lpmd/scalartable.h>
#include <lpmd/temporalproperty.h>
#include <lpmd/plugin.h>

class Vacf: public lpmd::ScalarTable, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Vacf(std::string args);
  ~Vacf();

  //Metodos propios de modulo vacf
  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(const std::vector<lpmd::SimulationCell> & simcell, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
};

#endif

