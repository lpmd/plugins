//
//
//

#ifndef __CRYSTALFCC_H__
#define __CRYSTALFCC_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class FCCGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   //
   FCCGenerator(std::string args);
   virtual ~FCCGenerator();

   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   void Generate(lpmd::SimulationCell & sc) const;

 private:
   int spc;
   long nx, ny, nz;
   double a;
};

#endif

