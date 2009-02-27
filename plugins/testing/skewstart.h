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
  //Metodos Generales
  SkewStartGenerator(std::string args); 
  virtual ~SkewStartGenerator();
  void ShowHelp() const;

  //Metodos propios del modulo skewstart
  void Generate(lpmd::SimulationCell & sc) const;

 private:
  long n;   // Number of atoms
  int spc;  // which species (atomic number)
};

#endif


