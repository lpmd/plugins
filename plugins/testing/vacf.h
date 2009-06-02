//
//
//

#ifndef __VACF_H__
#define __VACF_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>

class Vacf: public lpmd::StoredValue<lpmd::Matrix>, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Vacf(std::string args);
  ~Vacf();

  //Metodos propios de modulo vacf
  void Evaluate(lpmd::SimulationHistory & hist, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
  double dt;
};

#endif

