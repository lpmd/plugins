//
//
//

#ifndef __DISPVOL_H__
#define __DISPVOL_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>



class DispVol: public lpmd::StoredValue<lpmd::Matrix>, public lpmd::TemporalProperty, public lpmd::Module
{
 public:
  DispVol(std::string args);
  ~DispVol();

  std::string Keywords() const;

  void Evaluate(lpmd::ConfigurationSet & hist, lpmd::Potential & pot);

 private:
  lpmd::Matrix * m;
  long int delta_t;
};

#endif

