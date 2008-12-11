/*
 *
 *
 *
 */

#include <string>
#include <lpmd/simulationcell.h>
#include <lpmd/matrix.h>
#include <lpmd/potential.h>

extern "C" const char * PluginVersion();

lpmd::Matrix* gdr(SimulationCell & simcell,Potential & pot,long int nb,double rcut);
lpmd::Matrix* vacf(const std::vector<SimulationCell> & simcell, Potential & pot, double dt); 
