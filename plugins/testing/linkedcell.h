//
//
//

#ifndef __LINKEDCELL_MODULE_H
#define __LINKEDCELL_MODULE_H

#include <lpmd/simulationcell.h>
#include <lpmd/cellmanager.h>
#include <lpmd/plugin.h>

#include <exception>

class SimulationCellTooSmall: public std::exception 
{
 const char * what() const throw();
};

//
//
//

class AtomItem
{
 public:
   long i;
   AtomItem *prev, *next;
  
   AtomItem(long index): i(index) { };
};

class BasicSubCell
{
 public:
   Vector center;
   long index;
   AtomItem * atomhead, *atomtail;

   BasicSubCell(): atomhead(NULL), atomtail(NULL) { ClearAtoms(); }
   ~BasicSubCell();

   void ClearAtoms();
   void AddAtom(long index);
   AtomItem * GetAtomList() const { return atomhead; }
};

class NeighborSubCell
{
 public:
   BasicSubCell * cell;
   NeighborSubCell *prev, *next;
   Vector disp;

   NeighborSubCell(): cell(NULL), prev(NULL), next(NULL), disp(Vector(0.0, 0.0, 0.0)) { }
   NeighborSubCell(BasicSubCell & scell, const Vector & v): cell(&scell), prev(NULL), next(NULL), disp(v) { }
};

class SubCell: public BasicSubCell
{
 public:
  SubCell();
  ~SubCell();

  void Allocate(long s);
  void AddFirstHalfNeighbor(BasicSubCell & cell, const Vector & disp);
  void AddSecondHalfNeighbor(BasicSubCell & cell, const Vector & disp);

  long GetNumberOfFirstHalf() const { return ifirst; }
  long GetNumberOfSecondHalf() const { return isecond; }

  NeighborSubCell & GetFirstHalfNeighbor(long i) const { return neighbors[i]; }
  NeighborSubCell & GetSecondHalfNeighbor(long i) const { return neighbors[iwall+i]; }

 private:
  long ifirst, isecond, iwall;
  NeighborSubCell *neighbors;
};

//
//
//
class LinkedCellManager
{
 public:
   LinkedCellManager(SimulationCell & cell, long nx, long ny, long nz, double rcut);
   ~LinkedCellManager();

   void BuildSubCellList(double rcut);
   double GetCutoff() const { return rcut; }
   void FillCells();
   long SubCellIndex(long * n, Vector & trans) const;
   long SubCellIndex(long i, long j, long k, Vector & trans) const;
   SubCell * GetSubCellList() const;

   long NumberOfSubCells() const { return grid[0]*grid[1]*grid[2]; }

   SubCell & GetSubCellByAtom(long i) const;

   SubCell & operator[](long i) const { return subcells[i]; }

   int GridSize(int i) const { return grid[i]; }

  private:
   long grid[3];
   double rcut;
   SimulationCell & realcell;
   SubCell * subcells;
};


class LinkedCellCellManager: public lpmd::CellManager, public lpmd::Module
{
 public:

   LinkedCellCellManager(std::string args);
   ~LinkedCellCellManager();

   void Show(std::ostream & os) const;
   std::string Keywords() const;

   void Reset();
   void UpdateCell(SimulationCell & sc);
   void BuildNeighborList(SimulationCell & sc, long i, std::list<Neighbor> & nlist, bool full);

 private:
   double rcut;
   long nx, ny, nz;
   LinkedCellManager * lcm;
};


#endif


