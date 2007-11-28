//
//
//

#ifndef __EWALD2_H__
#define __EWALD2_H__

#include <complex>

#include <lpmd/simulationcell.h>
#include <lpmd/potential.h>
#include <lpmd/plugin.h>

class Ewald2: public lpmd::Potential , public lpmd::Module
{
 public:
  /// Varialbes principales.
  double alpha,ep;            //parametro de ewald y epsilon'
  bool status;                //estado de ewald, inicializando (true), ya inicializado (false)
  double SE;                  //autoenergia (no depende de coordenadas.)
  double rc,kc;               //cutoffs real y reciproco.
  std::vector<lpmd::Vector> kpoints;//kpoints para la suma reciproca.
  //Settings de las variables publicas.
  void ReSet(lpmd::SimulationCell &sc);

  /// Constructor por omisión.
  Ewald2(std::string args);
  /// Destructor.
  ~Ewald2(){};

  /// Implementa Potential::energy.
  double energy(lpmd::SimulationCell & sc);
  /// Implementa Potential::UpdateForces .
  void UpdateForces(lpmd::SimulationCell & sc);

  /*Energias*/
  /// Retorna Termino del Espacio Real de la Energia de la celda.
  double RealEnergy(lpmd::SimulationCell & sc) const;
  /// Retorna Termino del Espacio Reciproco de la Energia de la celda.
  double ReciprocalEnergy(lpmd::SimulationCell & sc) const;
  /// Retorna Auto-Energia de la Celda.
  double SelfEnergy(lpmd::SimulationCell & sc) const;
  /// Retorna Correccion Dipolar de la Energia de la celda.
  double DipoleCorrectionEnergy(lpmd::SimulationCell & sc) const;

  /*Fuerzas*/
  /// Retorna la constribucion a la fuerza sobre el atomo i hecha por el termino del espacio real.
  lpmd::Vector RealForce(lpmd::SimulationCell & sc,int i) const;
  /// Retorna la contribucion a la fuerza sobre el atomo i hecha por el termino del espacio reciproco.
  lpmd::Vector ReciprocalForce(lpmd::SimulationCell & sc,int i) const;
  /// Retorna la contribucion a la fuerza sobre el atomo i hecha por la correcion dipolar.
  lpmd::Vector DipoleCorrectionForce(lpmd::SimulationCell & sc,int i) const;

  /*Funciones Extras*/
  void SetParameter(std::string name);
  void Show() const;
  std::string Keywords() const;
};

#endif
