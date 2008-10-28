//
//
//

#ifndef __SUTTON_H__
#define __SUTTON_H__

#include <lpmd/metalpotential.h>
#include <lpmd/plugin.h>
#include <lpmd/vector.h>

class SuttonChen: public lpmd::MetalPotential , public lpmd::Module
{
 public:
  //Metodos Generales
  SuttonChen(std::string args);
  ~SuttonChen(){};
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos Propios de modulo suttonchen
  double pairEnergy(const double &r) const;
  double rhoij(const double &r) const;
  double F(const double &rhoi) const;
  lpmd::Vector PairForce(const lpmd::Vector &rij) const;
  lpmd::Vector ManyBodies(const lpmd::Vector &rij, const double &rhoi, const double &rhoj) const;
  double deltarhoi(const double & rhobar, const int & N) const;
  double deltaU1(const double & rhobar, const int & N) const;
  double deltaU2(const double & rhobar, const int & N, const double & rhoi) const;
  double VirialCorrection(const double &rhobar, const int &N, const double &rhoi) const;

 private:
  double e,a,n,m,c,rcut;
};

#endif
