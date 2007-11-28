//
//
//

#ifndef __LJ_POTENTIAL_H__
#define __LJ_POTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>

class LennardJones: public lpmd::PairPotential, public lpmd::Module
{
 public:

   // LJ potential parameters
   double sigma, epsilon, cutoff;

   // Constructor y Destructor
   LennardJones(std::string args); 
   ~LennardJones() { };

  double pairEnergy(const double & r) const;
  lpmd::Vector pairForce(const lpmd::Vector & r) const;

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;

};


#endif

