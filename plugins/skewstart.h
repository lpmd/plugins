//
//
//

#ifndef __SKEWSTARTMODULE_H__
#define __SKEWSTARTMODULE_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class SkewStartGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   //
   SkewStartGenerator(std::string args); 
   virtual ~SkewStartGenerator();

   void Generate(lpmd::SimulationCell & sc) const;

   // from Module
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

 private:
   long n;   // Number of atoms
   int spc;  // which species (atomic number)
};

#endif


