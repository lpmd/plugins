//
//
//

#ifndef __CNA_H__
#define __CNA_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

#include <map>

class IndexTrio
{
 public:
  unsigned int j, k, l;

  IndexTrio(unsigned int j0, unsigned int k0, unsigned int l0): j(j0), k(k0), l(l0) { }
};

inline bool operator<(const IndexTrio & t1, const IndexTrio & t2)
{
 if (t1.j < t2.j) return true;
 if (t1.j > t2.j) return false;
 if (t1.k < t2.k) return true;
 if (t1.k > t2.k) return false;
 if (t1.l < t2.l) return true;
 return false;
}

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
    int mode, spc1, spc2;
    double rcut;
    std::map<IndexTrio, int> refmap;
};

#endif

