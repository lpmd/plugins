//
//
//

#ifndef __MSD_H__
#define __MSD_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/temporalproperty.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class MSD: public Value<Matrix>, public TemporalProperty, public Module
{
 public:
  MSD(std::string args);
  ~MSD();

  const Matrix & CurrentValue() const { return *m; }
  void Evaluate(const SimulationHistory & hist, Potential & pot);

 private:
  Matrix * m;
};

#endif

