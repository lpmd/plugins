//
//
//

#ifndef __EWALD_H__
#define __EWALD_H__

#include <complex>

#include <lpmd/simulationcell.h>
#include <lpmd/potential.h>
#include <lpmd/plugin.h>

class Ewald: public lpmd::Potential , public lpmd::Module
{
 public:
  /// Varialbes principales.
  double alpha,ep;
  /// Constructor por omisi�n.
  Ewald(std::string args);
  /// Destructor.
  ~Ewald(){};

  /// Implementa Potential::energy .
  double energy(lpmd::SimulationCell & sc) ;
  /// Implementa Potential::UpdateForces .
  void UpdateForces(lpmd::SimulationCell & sc);

  /// Retorna Termino del Espacio Real de la Energia de la celda.
  double RealEnergy(lpmd::SimulationCell & sc) const;
  /// Retorna Termino del Espacio Reciproco de la Energia de la celda.
  double ReciprocalEnergy(lpmd::SimulationCell &sc) const;
  /// Retorna Auto-Energia de la Celda.
  double SelfEnergy(lpmd::SimulationCell &sc) const;
  /// Retorna Correccion Dipolar dela Energia de la celda.
  double DipoleCorrectionEnergy(lpmd::SimulationCell &sc) const;

  /// Retorna la constribucion a la fuerza sobre el atomo i hecha por el termino del espacio real.
  lpmd::Vector RealForce(lpmd::SimulationCell &sc,int i) const;
  /// Retorna la contribucion a la fuerza sobre el atomo i hecha por el termino del espacio reciproco.
  lpmd::Vector ReciprocalForce(lpmd::SimulationCell &sc,int i) const;
  /// Retorna la contribucion a la fuerza sobre el atomo i hecha por la correcion dipolar.
  lpmd::Vector DipoleCorrectionForce(lpmd::SimulationCell &sc,int i) const;

  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;
};

#endif
