//
//
//

#ifndef __DISPVOL_H__
#define __DISPVOL_H__

#include <lpmd/value.h>
#include <lpmd/matrix.h>
#include <lpmd/temporalproperty.h>
#include <lpmd/plugin.h>

class DispVol: public lpmd::Value<lpmd::Matrix>, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  DispVol(std::string args);
  ~DispVol();

  std::string Keywords() const;

  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(const std::vector<lpmd::SimulationCell> & simcell, lpmd::Potential & pot);

 private:
  std::string outputfile;
  lpmd::Matrix * m;
  long int delta_t;
};

#endif

