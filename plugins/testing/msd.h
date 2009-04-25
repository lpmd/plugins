//
//
//

#ifndef __MSD_H__
#define __MSD_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/temporalproperty.h>
#include <lpmd/plugin.h>

class MSD: public lpmd::Value<lpmd::Matrix>, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  MSD(std::string args);
  ~MSD();

  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(const std::vector<lpmd::SimulationCell> & simcell, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
};

#endif

