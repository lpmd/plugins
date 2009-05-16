//
//
//

#ifndef __CRYSTALFCC_H__
#define __CRYSTALFCC_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class FCCGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   //Metodos Generales
   FCCGenerator(std::string args);
   virtual ~FCCGenerator();
   void ShowHelp() const;

   //Metodos Propios del modulo fccgenerator
   void Generate(Configuration & config) const;

 private:
   int spc;
   long nx, ny, nz;
};

#endif
