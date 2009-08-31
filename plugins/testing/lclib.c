/*
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "lclib.h"

inline long int wrap(long int i, long int m)
{
 if (i < 0) return (i+m);
 if (i > m-1) return (i-m);
 return i;
}

void LC_Init(LinkedCell * lc, int nx, int ny, int nz)
{
 lc->grid[0] = nx;
 lc->grid[1] = ny;
 lc->grid[2] = nz;
 lc->cells = NULL;
 lc->dist2 = NULL;
 lc->atomvol = 0;
 for (int q=0;q<3;++q) lc->cdtable[q] = (double *)(malloc((lc->grid[q])*(2*(lc->grid[q])+1)*sizeof(double)));
}

void LC_FillCells(LinkedCell * lc, unsigned int natoms, double * pos, double * box)
{
 long int ncells = (lc->grid[0])*(lc->grid[1])*(lc->grid[2]);
 int maxatoms = (int)(3*((double)(natoms))/ncells)+1;
 if (maxatoms < 20) maxatoms = 20;
 if (lc->cells == NULL)
 {
  lc->cells = (InnerCell *)(malloc(ncells*sizeof(InnerCell)));
  for (int i=0;i<ncells;++i) 
  {
   (lc->cells[i]).atoms = (long int *)(malloc(maxatoms*sizeof(long int)));
   (lc->cells[i]).neighcells = NULL;
   (lc->cells[i]).ncda = NULL;
  }
 }
 int i, j, k;
 InnerCell * this_cell;
 for (i=0;i<ncells;++i) 
 {
  this_cell = &((lc->cells)[i]);
  for (int q=0;q<maxatoms;++q) this_cell->atoms[q] = -1;
  this_cell->n = 0;
 }
 unsigned long int n=0;
 for (n=0;n<natoms;++n)
 {
  i = (int)(floor((lc->grid[0])*pos[n*3]/box[0]));
  j = (int)(floor((lc->grid[1])*pos[n*3+1]/box[1]));
  k = (int)(floor((lc->grid[2])*pos[n*3+2]/box[2]));
  //fprintf(stderr, "DEBUG atom %d, <%8.8f %8.8f %8.8f>\n", n, pos[n*3], pos[n*3+1], pos[n*3+2]);
  /*
  assert(pos[n*3] >= 0.0);
  assert(pos[n*3] <= box[0]);
  assert(pos[n*3+1] >= 0.0);
  assert(pos[n*3+1] <= box[1]);
  assert(pos[n*3+2] >= 0.0);
  assert(pos[n*3+2] <= box[2]);
  */
  if (i < 0) i += lc->grid[0];
  else if (i > lc->grid[0]-1) i -= lc->grid[0];
  if (j < 0) j += lc->grid[1];
  else if (j > lc->grid[1]-1) j -= lc->grid[1];
  if (k < 0) k += lc->grid[2];
  else if (k > lc->grid[2]-1) k -= lc->grid[2];
  long int cind = k+j*(lc->grid[2])+i*(lc->grid[2])*(lc->grid[1]);
  //fprintf(stderr, "DEBUG cind=%d\n", cind);
  //assert((cind >= 0) && (cind <= ncells-1));
  this_cell = &((lc->cells)[cind]);
  //assert(this_cell->n < maxatoms-1);
  this_cell->atoms[this_cell->n] = n;
  this_cell->atoms[(this_cell->n)+1] = -1;
  (this_cell->n)++;
 }
 lc->atomvol = (box[0]*box[1]*box[2])/natoms;
 
 // Precalculates all the cell displacements
 for (int q=0;q<3;++q)
   for (int ik=0;ik<lc->grid[q];++ik)
      for (int dd=-(lc->grid[q]-1);dd<=(lc->grid[q]-1);++dd)
      {
       int p = ik*(2*lc->grid[q]+1)+(dd+lc->grid[1]-1);
       lc->cdtable[q][p] = 0.0;
       if (ik + dd >= lc->grid[q]-1) lc->cdtable[q][p] = box[q]*floor((ik+dd)/(lc->grid[q]-1));
       else lc->cdtable[q][p] = -box[q]*(floor((-ik-dd-1)/(lc->grid[q]-1))+1);
      }

 //printf("DEBUG maxatoms=%d\n", maxatoms);
 for (i=0;i<ncells;++i)
 {
  this_cell = &((lc->cells[i]));
  //printf("DEBUG cell %d has %lu atoms\n", i, this_cell->n);
  //printf("DEBUG terminator: %d\n", this_cell->atoms[this_cell->n]);
  //assert(this_cell->atoms[this_cell->n] == -1);
 }
}

void LC_SetDistanceSqrFunction(LinkedCell * lc, double (*distfunc)(double,double,double))
{
 lc->dist2 = distfunc; 
 //fprintf(stderr, "DEBUG calling LC_SetDistanceSqrFunction\n");
}

unsigned long int LC_NumberOfAtomsInCell(LinkedCell * lc, unsigned int ci)
{
 InnerCell * icell =&((lc->cells)[ci]);
 long int * atom; 
 unsigned long int s=0;
 for (atom = icell->atoms;(*atom) >= 0;atom++) s++;
 return s;
}

long int * LC_GetMINeighbors(LinkedCell * lc, unsigned long int n, double rcut, unsigned long int natoms, double * pos, double * box)
{
 assert(box > (void *)NULL);
 assert(&lc != 0); //icc 869
 int i, nn=0;
 int enn = (int)((4.0/3.0)*3.141592653589793*pow(2.0*rcut, 3.0)/(lc->atomvol));
 long int * nlist = (long int *)(malloc(2*enn*sizeof(long int)));
 for (i=0;i<natoms;++i)
 {
  if (i == n) continue;
  double rr2 = (*(lc->dist2))(pos[n*3]-pos[i*3], pos[n*3+1]-pos[i*3+1], pos[n*3+2]-pos[i*3+2]);
  if (rr2 < rcut*rcut) nlist[nn++] = i;
 }
 nlist[nn] = -1;
 return nlist;
}

inline void ProcessIntraAtom(int n, NeighborInfo * nlist, int * nn, long int * atom, double rcut, double * pos)
{
// fprintf(stderr, "DEBUG ProcessIntraAtom\n");
 int q;
 double dr[3];
 for (q=0;q<3;++q) dr[q] = pos[n*3+q]-pos[(*atom)*3+q];
 double rr2 = (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
 if (rr2 < rcut*rcut) 
 {
  //fprintf(stderr, "DEBUG intra*** rr=%8.8f\n", sqrt(rr2));
  nlist[*nn].i = n;
  nlist[*nn].j = (*atom);
  nlist[*nn].r2 = rr2;
  for (q=0;q<3;++q) nlist[*nn].dr[q] = dr[q];
  (*nn)++;
 }
}

inline void GetCellDisplacement(LinkedCell * lc, double * da, int * index, int * del)
{
 for (int q=0;q<3;++q)
     da[q] = lc->cdtable[q][index[q]*(2*lc->grid[q]+1)+(del[q]+lc->grid[1]-1)];
}

inline double Distance2(int n, double * dr, double * da, long int * atom, double * pos, double * box)
{
 double s = 0.0;
 for (int q=0;q<3;++q)
 {
  dr[q] = pos[n*3+q]-pos[(*atom)*3+q]-da[q];
  if (dr[q] >= 0.5*box[q]) dr[q] -= box[q];
  else if (dr[q] < -0.5*box[q]) dr[q] += box[q];
  s += dr[q]*dr[q];
 }
 return s;
}

inline int ProcessInterAtom(int n, NeighborInfo * nlist, LinkedCell * lc, int * nn, long int * atom, double * dr, double rr2, double rcut)
{
 assert(&lc != 0); // icc 869
 if ((n == *atom) && (rr2 < 1.0E-08)) return 0;
 else
 {
  if (rr2 < rcut*rcut)
  {
   //fprintf(stderr, "DEBUG inter!!! rr=%8.8f\n", sqrt(rr2));
   nlist[*nn].i = n;
   nlist[*nn].j = (*atom);
   nlist[*nn].r2 = rr2;
   for (int q=0;q<3;++q) nlist[*nn].dr[q] = (-dr[q]); // minus is needed for lpmd...
   (*nn)++;
   return 1;
  }
 }
 return 0;
}

inline void GenerateSubcellsList(LinkedCell * lc, int * index, double rcut, double * box)
{
 int maxdx = 1+(int)(floor((rcut/(box[0]/(lc->grid[0])))));
 int maxdy = 1+(int)(floor((rcut/(box[1]/(lc->grid[1])))));
 int maxdz = 1+(int)(floor((rcut/(box[2]/(lc->grid[2])))));

 if (lc->grid[0] < 3) maxdx = 1;
 if (lc->grid[1] < 3) maxdy = 1;
 if (lc->grid[2] < 3) maxdz = 1;

 if (maxdx < 0) maxdx = 0;
 if (maxdy < 0) maxdy = 0;
 if (maxdz < 0) maxdz = 0;
 
 double delta[3]; 
 for (int q=0;q<3;++q) delta[q] = box[q]/(lc->grid[q]); 

 double rcrit = rcut+sqrt(pow(delta[0], 2.0)+pow(delta[1], 2.0)+pow(delta[2], 2.0));

 InnerCell * this_cell = &((lc->cells)[index[2]+index[1]*(lc->grid[2])+index[0]*(lc->grid[2])*(lc->grid[1])]);
 //assert(this_cell->neighcells == NULL);
 this_cell->neighcells = (InnerCell **)(malloc((2*maxdx+1)*(2*maxdy+1)*(2*maxdz+1)*sizeof(InnerCell *)));
 this_cell->ncda = (double *)(malloc(3*(2*maxdx+1)*(2*maxdy+1)*(2*maxdz+1)*sizeof(double)));

 int nc=0, del[3];
 for (del[0]=-maxdx;del[0]<=maxdx;del[0]++)
   for (del[1]=-maxdy;del[1]<=maxdy;del[1]++)
     for (del[2]=-maxdz;del[2]<=maxdz;del[2]++)
     {
      double d2cell = pow(del[0]*delta[0],2.0)+pow(del[1]*delta[1],2.0)+pow(del[2]*delta[2],2.0);
      if (d2cell > rcrit*rcrit) continue;
      unsigned long int kkk = wrap(index[2]+del[2], lc->grid[2])+wrap(index[1]+del[1], lc->grid[1])*(lc->grid[2])+wrap(index[0]+del[0], lc->grid[0])*(lc->grid[2])*(lc->grid[1]);
      InnerCell * outer_cell = &((lc->cells)[kkk]);
      if (outer_cell == this_cell) continue;
      this_cell->neighcells[nc] = outer_cell; 
      double da[3];
      GetCellDisplacement(lc, da, index, del);
      for (int q=0;q<3;++q) this_cell->ncda[3*nc+q] = da[q];
      nc++;
     }
 this_cell->neighcells[nc] = NULL;
 //double cellvol = ((double)(box[0]*box[1]*box[2]))/(lc->grid[0]*lc->grid[1]*lc->grid[2]);
 //int prednc = (int)((4.0/3.0)*3.141592653589793*pow(rcrit, 3.0)/cellvol);
 //fprintf(stderr, "DEBUG nc = %d, predicted nc = %d\n", nc, prednc);
}

NeighborInfo * LC_GetNeighbors(LinkedCell * lc, unsigned long int n, double rcut, double * pos, double * box, int full)
{
// fprintf(stderr, "DEBUG calling LC_GetNeighbors\n");
 int index[3], nn=0;
 index[0] = (int)(floor((lc->grid[0])*pos[n*3]/box[0]));
 index[1] = (int)(floor((lc->grid[1])*pos[n*3+1]/box[1]));
 index[2] = (int)(floor((lc->grid[2])*pos[n*3+2]/box[2]));
 if (index[0] < 0) index[0] += lc->grid[0];
 else if (index[0] > lc->grid[0]-1) index[0] -= lc->grid[0];
 if (index[1] < 0) index[1] += lc->grid[1];
 else if (index[1] > lc->grid[1]-1) index[1] -= lc->grid[1];
 if (index[2] < 0) index[2] += lc->grid[2];
 else if (index[2] > lc->grid[2]-1) index[2] -= lc->grid[2];
 //fprintf(stderr, "DEBUG pos: %8.8f  %8.8f  %8.8f  (%8.8f %8.8f %8.8f)\n", pos[n*3], pos[n*3+1], pos[n*3+2], box[0], box[1], box[2]);
 /*
 assert(index[0] >= 0);
 assert(index[0] <= lc->grid[0]-1);
 assert(index[1] >= 0);
 assert(index[1] <= lc->grid[1]-1);
 assert(index[2] >= 0);
 assert(index[2] <= lc->grid[2]-1);
 */
 InnerCell * this_cell = &((lc->cells)[index[2]+index[1]*(lc->grid[2])+index[0]*(lc->grid[2])*(lc->grid[1])]);
 //fprintf(stderr, "DEBUG this_cell = %p\n", (void*)(this_cell));
 //fprintf(stderr, "DEBUG this_cell->neighcells = %p\n", (void*)(this_cell->neighcells));

 if (this_cell->neighcells == NULL)
 {
  GenerateSubcellsList(lc, index, rcut, box);
 }

 int enn = (int)((4.0/3.0)*(3.141592653589793)*pow(2.0*rcut, 3.0)/(lc->atomvol));
 NeighborInfo * nlist = (NeighborInfo *)(malloc(2*enn*sizeof(NeighborInfo)));

 // 
 // Fills the full or half list 
 //
 // Intracell first:
 long int * atom;
 //assert(this_cell != NULL);
 //assert(this_cell->atoms != NULL);
 for (atom = this_cell->atoms;(*atom) >= 0;atom++) 
 {
  if (*atom == n) continue;
  if (!full && (*atom > n)) continue;
  ProcessIntraAtom(n, nlist, &nn, atom, rcut, pos);
 }
 
 int nc = 0, st = -1, nwin=0, nfail=0;
 while (1)
 {
  InnerCell * outer_cell = this_cell->neighcells[nc];
  if (outer_cell == NULL) break;
  double * da = &(this_cell->ncda[3*nc]);
  double rr2, dr[3];
  for (atom = outer_cell->atoms;(*atom) >= 0;atom++)
  {
   if (!full && (*atom > n)) continue;
   rr2 = Distance2(n, dr, da, atom, pos, box);
   st = ProcessInterAtom(n, nlist, lc, &nn, atom, dr, rr2, rcut);
   if (st == 1) nwin++;
   else nfail++;
  }
  nc++;
 }
 
 //fprintf(stderr, "%d %d\n", (2*maxdx+1)*(2*maxdy+1)*(2*maxdz+1), nc);
 //fprintf(stderr, "DEBUG %ld %ld\n", nn-nn0, pia_cnt);
 //fprintf(stderr, "DEBUG fail %8.2f\n", 100.0*nfail/((double)(nwin+nfail)));
 nlist[nn].j = -1;
 return nlist;
}

void LC_Destroy(LinkedCell * lc)
{
 //fprintf(stderr, "DEBUG calling LC_Destroy\n");
 if (lc->cells != NULL)
 {
  int i, j, k;
  for (i=0;i<lc->grid[0];++i)
   for (j=0;j<lc->grid[1];++j)
     for (k=0;k<lc->grid[2];++k) 
     {
      InnerCell * this_cell = &((lc->cells)[k+j*(lc->grid[2])+i*(lc->grid[2])*(lc->grid[1])]);
      free(this_cell->atoms);
      free(this_cell->neighcells);
      free(this_cell->ncda);
     }
  free(lc->cells);
  lc->cells = NULL;
 }
 for (int q=0;q<3;++q) free(lc->cdtable[q]);
}

