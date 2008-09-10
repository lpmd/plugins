//
//
//

#ifndef __CNA_H__
#define __CNA_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class CommonNeighborAnalysis: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
  //Metodos Generales
  CommonNeighborAnalysis(std::string args);
  ~CommonNeighborAnalysis();
  void ShowHelp() const;
  std::string Keywords() const;

  const lpmd::Matrix & Value() const { return *m; }
  void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int mode;
    double rcut;
};

#endif

