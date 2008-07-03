//
//
//

#ifndef __CRYSTAL2D_H__
#define __CRYSTAL2D_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class Crystal2DGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
   
  //Metodos Generales
  Crystal2DGenerator(std::string args);
  virtual ~Crystal2DGenerator();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios del module scgnerator
  void Generate(lpmd::SimulationCell & sc) const;

 private:
   int spc;
   long nx, ny, nz;
   double a, b, gamma;
};

#endif

