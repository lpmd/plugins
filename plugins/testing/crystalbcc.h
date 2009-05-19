//
//
//

#ifndef __CRYSTALBCC_H__
#define __CRYSTALBCC_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class BCCGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   //Metodos Generales
   BCCGenerator(std::string args);
   virtual ~BCCGenerator();
   void ShowHelp() const;

   //Metodos Propios del modulo fccgenerator
   void Generate(Configuration & conf) const;

 private:
   int spc;
   long nx, ny, nz;
};

#endif
