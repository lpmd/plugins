//
//
//

#ifndef __MINIMUMIMAGE_MODULE_H
#define __MINIMUMIMAGE_MODULE_H

#include <lpmd/cellmanager.h>
#include <lpmd/plugin.h>

class MinimumImageCellManager: public lpmd::CellManager, public lpmd::Module
{
 public:

   MinimumImageCellManager(std::string args);
   ~MinimumImageCellManager();

   void Show(std::ostream & os) const;

   void Reset();
   void UpdateCell(SimulationCell & sc);
   void BuildNeighborList(SimulationCell & sc, long i, std::vector<Neighbor> & nlist, bool full, double);
   double Cutoff() const;

 private:
   double rcut;
};


#endif


