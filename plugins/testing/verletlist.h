//
//
//

#ifndef __VERLETLIST_MODULE_H
#define __VERLETLIST_MODULE_H

#include <lpmd/cellmanager.h>
#include <lpmd/plugin.h>

class VerletListCellManager: public lpmd::CellManager, public lpmd::Module
{
 public:

   VerletListCellManager(std::string args);
   ~VerletListCellManager();

   void Show(std::ostream & os) const;
   std::string Keywords() const;

   void Reset();
   void UpdateCell(SimulationCell & sc);
   double Cutoff() const;
   void BuildNeighborList(SimulationCell & sc, long i, std::list<Neighbor> & nlist, bool full, double a) {}
   void BuildList(SimulationCell & sc, bool full, double rcut);
//   void BuildNeighborList(SimulationCell & sc, long i, bool full, double);
   void PrintTimes() {}

 private:
   double cutoff;
   std::vector<lpmd::Vector> oldposition;
   std::list<lpmd::Neighbor> *neigh;
   bool evaluate;
   long calls;
};


#endif


