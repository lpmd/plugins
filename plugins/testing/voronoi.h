//
//
//

#ifndef __VORONOIMODULE_H__
#define __VORONOIMODULE_H__

#include <lpmd/cellgenerator.h>
#include <lpmd/plugin.h>

class VoronoiGenerator: public lpmd::CellGenerator, public lpmd::Module
{
 public:
  //Metodos Generales
  VoronoiGenerator(std::string args); 
  virtual ~VoronoiGenerator();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos propios del modulo voronoi
  void Generate(lpmd::SimulationCell & sc) const;

 private:
  long Ncell;			// Number of cells with atoms
  double a;				// Net parameter (size of the base cell)
  int spc;				// which species (atomic number)
  std::string type;	// type of base cell
};

#endif


