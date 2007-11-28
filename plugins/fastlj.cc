//
//
//

#include <iostream>

#include "fastlj.h"

using namespace lpmd;

FastLJ::FastLJ(std::string args): Module("fastlj") 
{ 
 bins = 500;
 ProcessArguments(args);
 Tabulate();
}

FastLJ::~FastLJ() 
{ 
 delete [] etable;
 delete [] fftable;
}

void FastLJ::Tabulate()
{
 double r, e, ff, r6, r12;
 etable = new double[bins];
 fftable = new double[bins];
 std::cerr << "-> Building table for Fast LJ potential, " << bins << " bins" << '\n';
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
 std::cerr << "-> Table finished." << '\n';
}

double FastLJ::pairEnergy(const double & r) const
{
 if (r >= cutoff) return 0.0;
 // Interpolate the value of the energy from the table
 long k = long(floor(bins*r/cutoff));
 return etable[k];
}

Vector FastLJ::pairForce(const Vector & r) const
{
 double rr2 = r.Mod2();
 if (rr2 >= cutoff*cutoff) return Vector(0.0);
 double rr = sqrt(rr2);
 Vector fv = r;
 long k = long(floor(bins*rr/cutoff));
 fv.Scale(fftable[k]);
 return fv;
}

//
//
//
void FastLJ::SetParameter(std::string name)
{
 if (name == "sigma") 
 {
  AssignParameter("sigma", GetNextWord());
  sigma = GetDouble("sigma");
 }
 if (name == "epsilon")
 {
  AssignParameter("epsilon", GetNextWord());
  epsilon = GetDouble("epsilon");
 }
 if (name == "cutoff") 
 {
  AssignParameter("cutoff", GetNextWord());
  cutoff = GetDouble("cutoff");
 }
 if (name == "bins")
 {
  AssignParameter("bins", GetNextWord());
  bins = GetInteger("bins");
 }
}

void FastLJ::Show() const
{
 Module::Show();
 std::cout << "     sigma = " << sigma << '\n';
 std::cout << "   epsilon = " << epsilon << '\n';
 std::cout << "    cutoff = " << cutoff << '\n';
 std::cout << "      bins = " << bins << '\n';
}

std::string FastLJ::Keywords() const { return "sigma epsilon cutoff bins"; }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new FastLJ(args); }
void destroy(Module * m) { delete m; }


