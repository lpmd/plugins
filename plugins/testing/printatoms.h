//
//
//

#ifndef __PRINTATOMS_SM_H__
#define __PRINTATOMS_SM_H__

#include <lpmd/visualizer.h>
#include <lpmd/plugin.h>

class PrintAtomsVisualizer: public lpmd::Visualizer, public lpmd::Module
{
 public:

  //Metodos Generales 
  PrintAtomsVisualizer(std::string args);
  ~PrintAtomsVisualizer();
  void ShowHelp() const;

  //Metodos Propios de modulo printatoms
  void Apply(const lpmd::Simulation & sim);

 private:
  long int from_at, to_at;
};

#endif



