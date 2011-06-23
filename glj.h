//
//
//

#ifndef __GLJ_POTENTIAL_H__
#define __GLJ_POTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>

class GenericLennardJones: public lpmd::PairPotential, public lpmd::Plugin
{
 public:
  //Metodos Generales
  GenericLennardJones(std::string args); 
  ~GenericLennardJones() { };
  void ShowHelp() const;

  double pairEnergy(const double & r) const;
  lpmd::Vector pairForce(const lpmd::Vector & r) const;

 private:
  double sigma, epsilon, a, b, n, m, cutoff;
};


#endif

