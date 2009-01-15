//
//
//

#ifndef __EWALDPOTENTIAL_H__
#define __EWALDPOTENTIAL_H__

#include <lpmd/potential.h>
#include <lpmd/plugin.h>

#include <lpmd/vector.h>

class Ewald: public lpmd::Potential, public lpmd::Module
{
 public:
  //Metodos Generales
  Ewald(std::string args); 
  ~Ewald();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Proopios de modulo ewald
  double energy(lpmd::SimulationCell & sc);
  void UpdateForces(lpmd::SimulationCell & sc);

 private:
   int kmax;
   bool surfdip;
   double alpha, etol, rcut, kcut;
   std::vector<lpmd::Vector> * kpoints;
   //
   void BuildKPointMesh(lpmd::SimulationCell & sc);
   void RealSpace(lpmd::SimulationCell & sc, lpmd::Vector * forces, double & e);
   void ReciprocalSpace(lpmd::SimulationCell & sc, lpmd::Vector * forces, double & e);
   void SurfaceDipole(lpmd::SimulationCell & sc, lpmd::Vector * forces, double & e);
   double EnergyConstantCorrection(lpmd::SimulationCell & sc);
};


#endif

