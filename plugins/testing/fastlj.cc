//
//
//

#include "fastlj.h"

#include <iostream>

using namespace lpmd;

FastLJ::FastLJ(std::string args): Module("fastlj") 
{ 
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("sigma");
 DefineKeyword("epsilon");
 DefineKeyword("cutoff");
 DefineKeyword("bins", "500");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 sigma = GetDouble("sigma");
 epsilon = GetDouble("epsilon");
 cutoff = GetDouble("cutoff");
 SetCutoff(cutoff);
 bins = GetInteger("bins");
 Tabulate();
}

FastLJ::~FastLJ() 
{ 
 delete [] etable;
 delete [] fftable;
}

void FastLJ::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo implementa el potencial de Lennard Jones para interaccion de   \n";
 std::cout << " de pares.                                                                     \n";
 std::cout << "      Se utiliza la pairpotential de la API para llevar a cabo el calculo,     \n";
 std::cout << " ademas de generar una tabla con los valores para rapido acceso al potencial.  \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      sigma         : Especifica el valor de sigma para el potencial.          \n";
 std::cout << "      epsilon       : Especifica el valor para epsilon del potencial.          \n";
 std::cout << "      cutoff        : Radio de corte para el potencial interatomico.           \n";
 std::cout << "      bins          : Numero de celdas para tabla de potencial.                \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use fastlj as FLJARAR                                                         \n";
 std::cout << "     sigma 3.4                                                                 \n";
 std::cout << "     epsilon 0.5                                                               \n";
 std::cout << "     bins 500                                                                  \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential FLJARAR Ar Ar                                                     \n\n";
 std::cout << "      De esta forma seteamos el potencial de Lennard Jones entre los atomos    \n";
 std::cout << " de Ar con las constantes usadas en FLJARAR.                                   \n";
}

void FastLJ::Tabulate()
{
 double r, e, ff, r6, r12;
 etable = new double[bins];
 fftable = new double[bins];
 DebugStream() << "-> Building table for Fast LJ potential, " << bins << " bins" << '\n';
 for (long i=0;i<bins;++i)
 {
  r = cutoff*double(i)/double(bins);
  r6 = pow(sigma / r, 6.0e0);
  r12 = r6*r6;
  e = 4.0e0*epsilon*(r12 - r6); 
  ff = -48.0e0*(epsilon/(r*r))*(r12 - 0.50e0*r6);
  etable[i] = e;
  fftable[i] = ff;
 }
 DebugStream() << "-> Table finished." << '\n';
}

double FastLJ::pairEnergy(const double & r) const
{
 // Interpolate the value of the energy from the table
 long k = long(floor(bins*r/cutoff));
 return etable[k];
}

Vector FastLJ::pairForce(const Vector & r) const
{
 double rr2 = r.SquareModule();
 double rr = sqrt(rr2);
 Vector fv = r;
 long k = long(floor(bins*rr/cutoff));
 fv = fv*fftable[k];
// fv.Scale(fftable[k]);
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new FastLJ(args); }
void destroy(Module * m) { delete m; }
