//
//
//

#ifndef __NULL_PAIRPOTENTIAL_H__
#define __NULL_PAIRPOTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>

class NullPairPotential: public lpmd::PairPotential, public lpmd::Module
{
 public:

   // Constructor y Destructor
   NullPairPotential(std::string args); 
   ~NullPairPotential() { };

  double pairEnergy(const double & r) const;
  lpmd::Vector pairForce(const lpmd::Vector & r) const;

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;

};


#endif

