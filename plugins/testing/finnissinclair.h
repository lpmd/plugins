//
//
//

#ifndef __FINNIS_H__
#define __FINNIS_H__

#include <lpmd/metalpotential.h>
#include <lpmd/plugin.h>
#include <lpmd/vector.h>

class FinnisSinclair: public lpmd::MetalPotential , public lpmd::Plugin
{
 public:
  //Metodos Generales
  FinnisSinclair(std::string args);
  ~FinnisSinclair(){};
  void ShowHelp() const;

  //Metodos Propios de modulo FinnisSinclair
  double pairEnergy(const double &r) const;
  double rhoij(const double &r) const;
  double deltarhoi(const double&) const { return 0.0; }
  double F(const double &rhoi) const;
  lpmd::Vector PairForce(const lpmd::Vector &modrij, const double &mod) const;
  lpmd::Vector ManyBodies(const lpmd::Vector &modrij, const double &invrhoi, const double &invrhoj, const double & mod) const;
  double deltarhoi(const double & rhobar, const int & N) const;
  double deltaU1(const double & rhobar, const int & N) const;

 private:
  double c0,c1,c2,A,B,c,d;	// Parameters "c" and "d" are cutoffs
};

#endif
