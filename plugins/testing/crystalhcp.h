//
//
//

#ifndef __CRYSTALHCP_H__
#define __CRYSTALHCP_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class HCPGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   //Metodos Generales
   HCPGenerator(std::string args);
   virtual ~HCPGenerator();
   void ShowHelp() const;
   std::string Keywords() const;

   //Metodos Propios del modulo crystalhcp
   void Generate(lpmd::SimulationCell & sc) const;

 private:
   int spc;
   long nx, ny, nz;
};

#endif
