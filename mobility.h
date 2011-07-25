//
//
//

#ifndef __MOBILITY_H__
#define __MOBILITY_H__

#include <lpmd/storedvalue.h>
#include <lpmd/property.h>
#include <lpmd/matrix.h>
#include <lpmd/plugin.h>

class Mobility: public lpmd::StoredValue<lpmd::Matrix>, public lpmd::TemporalProperty, public lpmd::Plugin
{
 public:
  Mobility(std::string args);
  ~Mobility();
  void ShowHelp() const;

  void Evaluate(lpmd::ConfigurationSet & hist, lpmd::Potential & pot);

 private:
  std::string outputfile;
  long int delta_t;
  double rcutmin, rcutmax;
};

#endif

