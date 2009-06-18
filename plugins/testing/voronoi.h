//
//
//

#ifndef __VORONOIMODULE_H__
#define __VORONOIMODULE_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class VoronoiGenerator: public lpmd::CellGenerator, public lpmd::Plugin
{
 public:
  //Metodos Generales
  VoronoiGenerator(std::string args); 
  virtual ~VoronoiGenerator();
  void ShowHelp() const;

  //Metodos propios del modulo voronoi
  void Generate(lpmd::Configuration & conf) const;

 private:
  int spc;              // which species (atomic number)
  std::string type;     // type of base cell
  double a;             // Lattice constant (size of the base cell of each grain)
  int grains;          // Number of cells (grains) in to put in the configuration
};

#endif


