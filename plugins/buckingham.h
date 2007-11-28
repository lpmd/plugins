//
//
//

#ifndef __BUCK_POTENTIAL_H__
#define __BUCK_POTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>
//#include <lpmd/coulomb.h>

class Buckingham: public lpmd::PairPotential , public lpmd::Module
{
 public:

// Buckingham potential parameters
  double B1, Ro, B2, cutoff;
// Constructor y Destructor
   Buckingham(std::string args); 
   ~Buckingham() { };

  double pairEnergy(const double & r) const;
  lpmd::Vector pairForce(const lpmd::Vector & r) const;

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;
};


#endif

