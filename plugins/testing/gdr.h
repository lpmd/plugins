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
  void Evaluate(lpmd::Configuration & conf, lpmd::Potential & pot);

 private:
    double rcut;
    int nb;
    bool do_average;
};

#endif

