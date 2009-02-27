//
//
//

#ifndef __CRYSTALSC_H__
#define __CRYSTALSC_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class SCGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   
  //Metodos Generales
  SCGenerator(std::string args);
  virtual ~SCGenerator();
  void ShowHelp() const;

  //Metodos Propios del module scgnerator
  void Generate(lpmd::SimulationCell & sc) const;

 private:
   int spc;
   long nx, ny, nz;
};

#endif

