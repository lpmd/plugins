//
//
//

#ifndef __FASTLJ_POTENTIAL_H__
#define __FASTLJ_POTENTIAL_H__

#include <lpmd/pairpotential.h>
#include <lpmd/plugin.h>

class FastLJ: public lpmd::PairPotential, public lpmd::Module
{
 public:
   // Constructor y Destructor
   FastLJ(std::string args); 
   ~FastLJ();

   double pairEnergy(const double & r) const;
   lpmd::Vector pairForce(const lpmd::Vector & r) const;

   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

 private:
   long bins;
   double sigma, epsilon, cutoff;
   double *etable, *fftable;
   
   void Tabulate();
};


#endif

