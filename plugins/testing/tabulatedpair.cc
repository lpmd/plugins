//
//
//

#include "tabulatedpair.h"

#include <iostream>
#include <fstream>

using namespace lpmd;

TabulatedPair::TabulatedPair(std::string args): Plugin("fastlj", "2.0")
{ 
 ParamList & params = (*this);
 //
 DefineKeyword("file");
 DefineKeyword("cutoff");
 DefineKeyword("bins","500");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 cutoff = params["cutoff"];
 SetCutoff(cutoff);
 bins = int(params["bins"]);
 file = params["file"];
 ReadTable();
}

TabulatedPair::~TabulatedPair() 
{ 
 delete [] etable;
 delete [] ftable;
}

void TabulatedPair::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo lee un potencial tabulado de pares a partir de un fichero      \n";
 std::cout << " de datos.                                                                     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      cutoff         : Especifica el valor de radio de corte para el potencial.\n";
 std::cout << "      file           : Especifica el nombre del archivo.                       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use tabulatedpair as TABLJ                                                    \n";
 std::cout << "     file lennardjones.dat                                                     \n";
 std::cout << "     cutoff 1.90                                                               \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " potential TABLJ Ar Ar                                                       \n\n";
 std::cout << "      De esta forma seteamos el potencial de Lennard Jones tabulado entre los  \n";
 std::cout << " atomos de Ar.                                                                 \n";
}

void TabulatedPair::ReadTable()
{
 DebugStream() << "-> Opening Table-Potential filetype '"<<file<<"'."<<'\n';
 std::ifstream input;
 input.open(file.c_str());
 if(!input.is_open())
 {
  EndWithError("Impossible open the potential data file.");
 }

 //count the number lines ignore #.
 int n = 0;
 std::string line;
 while(!input.eof()) 
 {
  getline(input,line);
  RemoveUnnecessarySpaces(line);
  lpmd::Array<std::string> words = StringSplit(line,' ');
  DebugStream() << " -> line = " << line <<" size = " << words.Size() << '\n';
  if(line[0]!='#' && words.Size()==2) n++;
 }
 input.clear();input.seekg(0);
 DebugStream() << "-> Readed "<<n<<" lines " << '\n';
 bins = n;
 etable = new double[bins];
 ftable = new double[bins];
 int i=0;
 //Asign values
 while(!input.eof())
 {
  getline(input,line);
  RemoveUnnecessarySpaces(line);
  lpmd::Array<std::string> words = StringSplit(line,' ');
  if(words.Size()==2 && line[0]!='#')
  {
   etable[i] = atof(words[0].c_str());
   ftable[i] = atof(words[1].c_str());
   i++;
  }
 }
 DebugStream() << "-> Table was readed correctly." << '\n';
}

double TabulatedPair::pairEnergy(const double & r) const
{
 // Interpolate the value of the energy from the table
 long k = long(floor(bins*r/cutoff));
 return etable[k];
}

Vector TabulatedPair::pairForce(const Vector & r) const
{
 double rr2 = r.SquareModule();
 double rr = sqrt(rr2);
 Vector fv = r;
 long k = long(floor(bins*rr/cutoff));
 fv = fv*ftable[k];
 return fv;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new TabulatedPair(args); }
void destroy(Plugin * m) { delete m; }
