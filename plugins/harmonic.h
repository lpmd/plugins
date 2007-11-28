//
//
//

#ifndef __HARMONIC_POTENTIAL_H__
#define __HARMONIC_POTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>

class Harmonic: public lpmd::PairPotential, public lpmd::Module
{
 public:
   // Harmonic potential parameters
   double k;
   double a;

   // Constructor y Destructor
   Harmonic(std::string args); 
   ~Harmonic() { };

  double pairEnergy(const double & r) const;
  lpmd::Vector pairForce(const lpmd::Vector & r) const;

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;

};


#endif

