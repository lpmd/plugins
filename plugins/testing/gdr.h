//
//
//

#ifndef __GDR_H__
#define __GDR_H__

#include <lpmd/storedvalue.h>
#include <lpmd/matrix.h>
#include <lpmd/property.h>
#include <lpmd/plugin.h>

class Gdr: public lpmd::StoredValue<lpmd::Matrix>, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  Gdr(std::string args);
  ~Gdr();
  void ShowHelp() const;

  //Metodos Propios de modulo gdr
  const lpmd::Matrix & CurrentValue() const { return *m; }
  void Evaluate(lpmd::Configuration & conf, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    double rcut;
    int nb;
    bool do_average;
};

#endif

