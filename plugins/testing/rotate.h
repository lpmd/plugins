//
//
//

#ifndef __ROTATE_SM_H__
#define __ROTATE_SM_H__

#include <lpmd/systemmodifier.h>
#include <lpmd/vector.h>
#include <lpmd/plugin.h>
#include <lpmd/configuration.h>

using namespace lpmd;

class RotateModifier: public lpmd::SystemModifier, public lpmd::Module
{
 public:
  //Metodos Generales
  RotateModifier(std::string args);
  ~RotateModifier();

  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios
//  void Apply(lpmd::SimulationCell & sc);
//  void Apply(lpmd::MD & md);
  void Apply(Configuration & conf);
  void Apply(Simulation & md);
 
 private:
  lpmd::Vector axis;
  double rotmat[3][3];                      // arbitrary rotation matrix
  double angle;
};

#endif

