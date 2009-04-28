/*
 *
 *
 *
 */

#include <string>
#include <lpmd/simulationcell.h>
#include <lpmd/matrix.h>
#include <lpmd/potential.h>
#include <lpmd/basicvector.h>
#include <lpmd/vector.h>

#ifndef _AZAR_UTILS_H_
#define _AZAR_UTILS_H_

#include <iostream>
#include <cstdlib>  
#include <cmath>
#include <ctime>
#include <unistd.h>

void randomize();
void randomize(unsigned int) ;
double dazar(double=0.0e0, double = 1.0e0);
int iazar(int, int) ;

#endif

extern "C" const char * PluginVersion();

lpmd::Matrix* gdr(SimulationCell & simcell,Potential & pot,long int nb,double rcut);
lpmd::Matrix* vacf(const std::vector<SimulationCell> & simcell, Potential & pot, double dt);

// Replicate replica la celda de simulacion SC dentro de la celda sc,
// ReplicateBase toma la celda replicada y la pone dentro de una mas grande
void Replicate(SimulationCell & sc, unsigned long nx, unsigned long ny, unsigned long nz);void Rotate(SimulationCell & sc);
void Rotate(SimulationCell & sc, lpmd::Vector rotate);
void ReplicateRotate(const SimulationCell basecell, lpmd::Vector &center, lpmd::Vector &CellColor, unsigned long na, unsigned long nb, unsigned long nc, lpmd::Vector rotate, SimulationCell & simcell);
