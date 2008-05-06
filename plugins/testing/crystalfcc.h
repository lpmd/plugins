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
   //Metodos Generales
   FCCGenerator(std::string args);
   virtual ~FCCGenerator();
   void ShowHelp() const;
   std::string Keywords() const;

   //Metodos Propios del modulo fccgenerator
   void Generate(lpmd::SimulationCell & sc) const;

 private:
   int spc;
   long nx, ny, nz;
};

#endif
