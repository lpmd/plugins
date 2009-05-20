/*
 *
 *
 *
 */

#ifndef __INNER_LC_H__
#define __INNER_LC_H__

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

void LC_Init(LinkedCell * lc, int nx, int ny, int nz);

void LC_SetDistanceSqrFunction(LinkedCell * lc, double (*distfunc)(double,double,double));

void LC_FillCells(LinkedCell * lc, unsigned int natoms, double * pos, double * box);

unsigned long int LC_NumberOfAtomsInCell(LinkedCell * lc, unsigned int ci);

NeighborInfo * LC_GetNeighbors(LinkedCell * lc, unsigned long int i, double rcut, double * pos, double * box, int full);

long int * LC_GetMINeighbors(LinkedCell * lc, unsigned long int i, double rcut, unsigned long int natoms, double * pos, double * box);

void LC_Destroy(LinkedCell * lc);

#endif


