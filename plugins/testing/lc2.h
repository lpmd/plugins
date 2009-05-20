//
//
//

#ifndef __LINKEDCELL2_MODULE_H
#define __LINKEDCELL2_MODULE_H

//
//
struct t_InnerCell
{
 long int * atoms; 
 unsigned long int n;
 struct t_InnerCell ** neighcells;
 double * ncda;
};
typedef struct t_InnerCell InnerCell;

//
//
typedef struct 
{
 double atomvol;
 int grid[3];
 InnerCell * cells;
 double (*dist2)(double dx, double dy, double dz);
 double * cdtable[3];
} LinkedCell;

typedef struct
{
 long int i;
 long int j;
 double dr[3];
 double r2;
} NeighborInfo;

//
//
//

extern "C" void LC_Init(LinkedCell * lc, int nx, int ny, int nz);

extern "C" void LC_SetDistanceSqrFunction(LinkedCell * lc, double (*distfunc)(double,double,double));

extern "C" void LC_FillCells(LinkedCell * lc, unsigned int natoms, double * pos, double * box);

extern "C" unsigned long int LC_NumberOfAtomsInCell(LinkedCell * lc, unsigned int ci);

extern "C" NeighborInfo * LC_GetNeighbors(LinkedCell * lc, unsigned long int i, double rcut, double * pos, double * box, int full);

extern "C" long int * LC_GetMINeighbors(LinkedCell * lc, unsigned long int i, double rcut, unsigned long int natoms, double * pos, double * box);

extern "C" void LC_Destroy(LinkedCell * lc);

#include <lpmd/configuration.h>
#include <lpmd/cellmanager.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class LinkedCellCellManager2: public lpmd::CellManager, public lpmd::Module
{
 public:

   LinkedCellCellManager2(std::string args);
   ~LinkedCellCellManager2();

   void Show(std::ostream & os) const;
   std::string Keywords() const { return "cutoff nx ny nz"; }

   void Reset();
   void UpdateCell(Configuration & sc);
   double Cutoff() const { return cutoff; }

   void BuildNeighborList(Configuration & sc, long i, NeighborList & nlist, bool full, double rcut);

  private:
    static double DistanceFunction(double dx, double dy, double dz);

    double * pos;
    double cutoff;
    int nx, ny, nz;
    LinkedCell lc;
};

#endif


