//
//
//

#include "nullmetalpotential.h"

#include <iostream>

using namespace lpmd;

NullMetalPotential::NullMetalPotential(std::string args): Plugin("nullmetalpotential", "2.0")
{ 
 //
 ProcessArguments(args); 
}

void NullMetalPotential::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa un potencial nulo de tipo embedded atom, para       \n";
 std::cout << " propositos de depuracion y benchmarking.                                      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use nullmetalpotential                                                        \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential nullmetalpotential Cu Cu                                            \n\n";
}

double NullMetalPotential::pairEnergy(const double &r) const { return 0.0; }

double NullMetalPotential::rhoij(const double &r) const { return 0.0; }

double NullMetalPotential::F(const double &rhoi) const { return 0.0; }

Vector NullMetalPotential::PairForce(const Vector &rij) const { return Vector(0.0, 0.0, 0.0); }

Vector NullMetalPotential::ManyBodies(const Vector &rij, const double &rhoi, const double &rhoj) const { return Vector(0.0, 0.0, 0.0); }

double NullMetalPotential::deltarhoi(const double &rhobar) const { return 0.0; }

double NullMetalPotential::deltaU1(const double &rhobar, const int &N) const { return 0.0; }

double NullMetalPotential::deltaU2(const double &rhobar, const int &N, const double &rhoi) const { return 0.0; }

double NullMetalPotential::VirialCorrection(const double &rhobar, const int &N, const double &rhoi) const { return 0.0; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) {return new NullMetalPotential(args);}
void destroy(Plugin * m) { delete m; }

