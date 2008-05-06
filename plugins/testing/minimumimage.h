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
   std::string Keywords() const;

   void Reset();
   void UpdateCell(SimulationCell & sc);
   void BuildNeighborList(SimulationCell & sc, long i, std::list<Neighbor> & nlist, bool full);

 private:
   double rcut;
};


#endif


