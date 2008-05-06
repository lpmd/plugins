//
//
//

#ifndef __MSD_H__
#define __MSD_H__

#include <lpmd/scalartable.h>
#include <lpmd/temporalproperty.h>
#include <lpmd/plugin.h>

class MSD: public lpmd::ScalarTable, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  MSD(std::string args);
  ~MSD();

  std::string Keywords() const;

  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(const std::vector<lpmd::SimulationCell> & simcell, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
};

#endif

