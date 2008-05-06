//
//
//

#ifndef __PRINTATOMS_H__
#define __PRINTATOMS_CM_H__

#include <lpmd/visualizer.h>
#include <lpmd/plugin.h>

class PrintAtomsVisualizer: public lpmd::Visualizer, public lpmd::Module
{
 public:

  //Metodos Generales 
  PrintAtomsVisualizer(std::string args);
  ~PrintAtomsVisualizer();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo printatoms
  void Apply(const lpmd::MD & md);

 private:
   long from_at, to_at;
};

#endif



