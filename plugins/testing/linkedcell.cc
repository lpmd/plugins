//
//
//

#include "linkedcell.h"
#include <lpmd/error.h>
#include <lpmd/configuration.h>

using namespace lpmd;

//
//
//

void BasicSubCell::ClearAtoms() 
{
 for (AtomItem * ak = GetAtomList();ak!=NULL;)
 {
  AtomItem * p = ak->next;
  delete ak;
  ak = p;
 }
 atomhead = atomtail = NULL; 
}

SubCell::SubCell(): BasicSubCell(), ifirst(0), isecond(0), iwall(0), neighbors(NULL) { }

SubCell::~SubCell() { delete [] neighbors; }

void SubCell::Allocate(long s)
{
 iwall = long(floor(s/2))+2;
 neighbors = new NeighborSubCell[2*iwall];
}

void SubCell::AddFirstHalfNeighbor(BasicSubCell & cell, const Vector & disp)
{
 //assert(ifirst < iwall);
 neighbors[ifirst] = NeighborSubCell(cell, disp);
 ifirst++;
}

void SubCell::AddSecondHalfNeighbor(BasicSubCell & cell, const Vector & disp)
{
 //assert((iwall+isecond) < 2*iwall); 
 neighbors[iwall+isecond] = NeighborSubCell(cell, disp);
 isecond++;
}

void BasicSubCell::AddAtom(long index)
{
 AtomItem * atomp = new AtomItem(index);
 atomp->prev = atomtail;
 atomp->next = NULL;
 if (atomtail != NULL) atomtail->next = atomp;
 else atomhead = atomp;
 atomtail = atomp; 
}

BasicSubCell::~BasicSubCell() { ClearAtoms(); }

LinkedCellManager::LinkedCellManager(Configuration & cell, long nx, long ny, long nz, double rcut): realcell(cell)
{
 grid[0] = nx;
 grid[1] = ny;
 grid[2] = nz;
 BuildSubCellList(rcut);
}

void LinkedCellManager::BuildSubCellList(double rc)
{
 rcut = rc;
 const long nsubc = NumberOfSubCells();
 if (nsubc == 0) throw PluginError("linkedcell", "The simulation box is too small");
 subcells = new SubCell[nsubc];
 BasicCell & cell = realcell.Cell();
 double ll =0.0;
 for (int j=0;j<3;++j)
   if (cell[j].Module()/double(GridSize(j)) > ll) ll = cell[j].Module()/double(GridSize(j));
 const double cell_cutoff = rcut + sqrt(3.0)*ll;
 long ncx = long(GridSize(0)*cell_cutoff/cell[0].Module());
 long ncy = long(GridSize(1)*cell_cutoff/cell[1].Module());
 long ncz = long(GridSize(2)*cell_cutoff/cell[2].Module());
 long nestim = long((4.0*M_PI/3.0)*NumberOfSubCells()*(pow(cell_cutoff+ll, 3.0))/cell.Volume());
 for (long i=0;i<nsubc;++i) subcells[i].Allocate(nestim+2);

 for (long k=0;k<GridSize(2);++k)
   for (long j=0;j<GridSize(1);++j)
     for (long i=0;i<GridSize(0);++i)
     {
      Vector tmp;
      long ni = SubCellIndex(i, j, k, tmp);
      SubCell & this_cell = subcells[ni];
      this_cell.index = ni;
      Vector cell_center = cell.Cartesian(Vector((i+0.5)/double(GridSize(0)), (j+0.5)/double(GridSize(1)), (k+0.5)/double(GridSize(2))));
      this_cell.center = cell_center;
      // 
      Vector trans;
      int smart_counter = 0;
      int nv = 0;

      for (int rr=-ncx;rr<=ncx;++rr)
       for (int qq=-ncy;qq<=ncy;++qq)
	    for (int pp=-ncz;pp<=ncz;++pp)
	    {
	     if (((pp == 0) && (qq == 0)) && (rr == 0)) continue;
         const Vector realdisp = cell.ScaleByCell(Vector(pp/double(GridSize(0)), qq/double(GridSize(1)), rr/double(GridSize(2))));
	     if (realdisp.Module() <= cell_cutoff)
	     {
          long nj = SubCellIndex(i+pp, j+qq, k+rr, trans);
	      trans = cell.ScaleByCell(trans);
	      nv++;
          if ((smart_counter % 2) == 0) this_cell.AddFirstHalfNeighbor(subcells[nj], trans);
	      else this_cell.AddSecondHalfNeighbor(subcells[nj], trans);
	     }
         smart_counter++;
	    }
       }
}

long LinkedCellManager::SubCellIndex(long * n, Vector & trans) const
{
 trans[0] = 0.0e0;
 trans[1] = 0.0e0;
 trans[2] = 0.0e0;
 for (int j=0;j<3;++j)
 {
  if (n[j] < 0)
  {
   while (n[j] < 0)
   {
    n[j] += GridSize(j);
    trans[j] = trans[j]-1.0;
   }
  }
  else if (n[j] >= GridSize(j)) 
  {
   while (n[j] >= GridSize(j))
   {
    n[j] -= GridSize(j);
    trans[j] = trans[j]+1.0;
   }
  }
 }
 return (n[2]*GridSize(0)*GridSize(1) + n[1]*GridSize(0) + n[0]);
}

long LinkedCellManager::SubCellIndex(long i, long j, long k, Vector & trans) const
{
 long n[3] = {i, j, k};
 long sci = SubCellIndex(n, trans);
 return sci;
}

LinkedCellManager::~LinkedCellManager() 
{
 if (subcells != NULL) delete [] subcells;
 subcells = NULL;
}

void LinkedCellManager::FillCells() 
{
 for (long i=0;i<NumberOfSubCells();++i) subcells[i].ClearAtoms();
 long int n = realcell.Atoms().Size();
 for (long int i=0;i<n;++i) GetSubCellByAtom(i).AddAtom(i);
 //long int totalcnt = 0;
 //std::cerr << "DEBUG Total subcells: " << NumberOfSubCells() << '\n';
 /*
 for (long int i=0;i<NumberOfSubCells();++i)
 {
  SubCell & subc = (*this)[i];
  int cnt = 0;
  for (AtomItem * ak = subc.GetAtomList();ak!=NULL;ak=ak->next) cnt++;
  // std::cerr << "DEBUG subcell " << i << " has " << cnt << " atoms \n";
  totalcnt += cnt;
 }
 */
 //assert(totalcnt == n);
}

SubCell * LinkedCellManager::GetSubCellList() const { return subcells; }

SubCell & LinkedCellManager::GetSubCellByAtom(long i) const
{
 BasicCell & cell = realcell.Cell();
 BasicParticleSet & atoms = realcell.Atoms();
 Vector tmp, r = cell.Fractional(atoms[i].Position());
 long k[3];
 for (int j=0;j<3;++j)
 {
  double d = 1.0 / double(GridSize(j));
  k[j] = long(floor(r[j]/d));
  if (k[j] < 0) k[j] = 0;
  else if (k[j] > (GridSize(j)-1)) k[j] = GridSize(j)-1;
 }
 long ci = SubCellIndex(k, tmp);
 return subcells[ci];
}

//
//
//

LinkedCellCellManager::LinkedCellCellManager(std::string args): Plugin("linkedcell", "2.0"), lcm(NULL) 
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("cutoff");
 DefineKeyword("nx", "7");
 DefineKeyword("ny", "7");
 DefineKeyword("nz", "7");
 // hasta aqui los parametros por omision
 ProcessArguments(args); 
 rcut = params["cutoff"];
 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);
 successperc = 0.0;
 nsuccess = 0;
}

LinkedCellCellManager::~LinkedCellCellManager() 
{ 
 DebugStream() << " *** LinkedCell Efficiency percentage: " << 100.0*successperc/double(nsuccess) << " % " <<'\n';
 if (lcm != NULL) delete lcm;
}

void LinkedCellCellManager::Show(std::ostream & os) const
{
 Module::Show(os);
 if (lcm != NULL) 
 {
  os << "   Simulation cell was divided in a " << lcm->GridSize(0) << "x" << lcm->GridSize(1) << "x" << lcm->GridSize(2) << " grid";
  os << " = " << lcm->NumberOfSubCells() << " subcells" << '\n';
 }
}

void LinkedCellCellManager::Reset()
{
 //if (lcm != NULL) delete lcm;
 //lcm = NULL;
}

void LinkedCellCellManager::UpdateCell(Configuration & sc) 
{ 
 if (lcm == NULL) lcm = new LinkedCellManager(sc, nx, ny, nz, rcut);
 lcm->FillCells(); 
}

double LinkedCellCellManager::Cutoff() const { return lcm->GetCutoff(); }

void LinkedCellCellManager::BuildNeighborList(Configuration & sc, long i, NeighborList & nlist, bool full, double rcut)
{
 BasicParticleSet & atoms = sc.Atoms();
 nlist.Clear();
 //lcm->FillCells();
 //const double rcc = lcm->GetCutoff();
 BasicAtom & this_atom = atoms[i];
 SubCell & subcell = lcm->GetSubCellByAtom(i);

 long int nwin = 0, nfail = 0;

 for (AtomItem * ak = subcell.GetAtomList();ak!=NULL;ak=ak->next)
 {
  if ( (full && (ak->i != i)) || ((! full) && (ak->i > i)) )
  {
   AtomPair nn;
   nn.i = &this_atom;
   nn.j = &atoms[ak->i];
   nn.rij = nn.j->Position() - this_atom.Position();
   nn.r = nn.rij.Module();
   if (nn.r < rcut)
   {
    nlist.Append(nn);
    nwin++;
   }
   else if (nn.r > rcut) nfail++;
  }
 }
 //
 for (long p=0;p<subcell.GetNumberOfFirstHalf();++p)
 {
  NeighborSubCell & nscell = subcell.GetFirstHalfNeighbor(p);
  BasicSubCell & scell = *(nscell.cell);
  const Vector disp = nscell.disp;
  for (AtomItem * ak = scell.GetAtomList();ak!=NULL;ak=ak->next)
  {
   AtomPair nn;
   nn.i = &this_atom;
   nn.j = &atoms[ak->i];
   const Vector newpos = nn.j->Position()+disp;
   nn.rij = newpos-this_atom.Position();
   nn.r = nn.rij.Module();
   if (nn.r < rcut)
   {
    nlist.Append(nn);
    nwin++;
   }
   else if (nn.r > rcut) nfail++;
  }
 }
 // 
 if (full)
 {
  for (long p=0;p<subcell.GetNumberOfSecondHalf();++p)
  {
   NeighborSubCell & nscell = subcell.GetSecondHalfNeighbor(p);
   BasicSubCell & scell = *(nscell.cell);
   const Vector disp = nscell.disp;
   for (AtomItem * ak = scell.GetAtomList();ak!=NULL;ak=ak->next)
   {
    AtomPair nn;
    nn.i = &this_atom;
    nn.j = &atoms[ak->i];
    const Vector newpos = nn.j->Position()+disp;
    nn.rij = newpos-this_atom.Position();
    nn.r = nn.rij.Module();
    if (nn.r < rcut)
    {
     nlist.Append(nn); 
     nwin++;
    }
    else if (nn.r > rcut) nfail++;
   }
  }
 }
 if ((nwin+nfail) > 0)
 {
  successperc += (double(nwin)/double(nwin+nfail));
  nsuccess++;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LinkedCellCellManager(args); }
void destroy(Plugin * m) { delete m; }

