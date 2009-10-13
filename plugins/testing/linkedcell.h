//
//
//

#ifndef __LINKEDCELL_MODULE_H
#define __LINKEDCELL_MODULE_H

#include <lpmd/configuration.h>
#include <lpmd/cellmanager.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class LinkedCell: public lpmd::CellManager, public lpmd::Plugin
{
 public:
   LinkedCell(std::string args);
   ~LinkedCell();

   void Reset();
   void UpdateCell(Configuration & conf);
   double Cutoff() const { return cutoff; }

   void BuildNeighborList(Configuration & conf, long i, NeighborList & nlist, bool full, double rcut);

  private:
    double cutoff;
    int nx, ny, nz;
    int * head, * tail, * subcell;
    long cells_inside;
    long * atomlist, *indexc;
    bool mode;
};

#endif

