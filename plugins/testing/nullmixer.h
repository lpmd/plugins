//
//
//

#ifndef __NULLMIXER_H__
#define __NULLMIXER_H__

#include <lpmd/systemmixer.h>
#include <lpmd/simulationcell.h>
#include <lpmd/plugin.h>

class NullMixer: public lpmd::SystemMixer, public lpmd::Module
{
 public:

  //Metodos Generales 
  NullMixer(std::string args);
  ~NullMixer();
  void ShowHelp() const;

  //Metodos Propios de modulo nullmixer
  SimulationCell Apply(lpmd::SimulationCell & sc1, lpmd::SimulationCell & sc2);
};

#endif

