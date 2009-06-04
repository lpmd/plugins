//
//
//

#ifndef __MSD_H__
#define __MSD_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class MSD: public StoredValue<Matrix>, public TemporalProperty, public Module
{
 public:
  MSD(std::string args);

  void Evaluate(ConfigurationSet & hist, Potential & pot);
};

#endif

