//
//
//

#ifndef __MORSE_POTENTIAL_H__
#define __MORSE_POTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>

class Morse: public lpmd::PairPotential, public lpmd::Module
{
 public:

   // Morse potential parameters
   double depth, a, re;

   // Constructor y Destructor
   Morse(std::string args); 
   ~Morse() { };

  double pairEnergy(const double & r) const;
  lpmd::Vector pairForce(const lpmd::Vector & r) const;

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;

};


#endif

