//
//
//

#include "lpmd.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>

using namespace lpmd;

LPMDFormat::LPMDFormat(std::string args): Module("lpmd")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 linecounter = new long int;
 DefineKeyword("file");
 DefineKeyword("each", "1");
 DefineKeyword("level", "0");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 rcell = GetBool("replacecell");
}

LPMDFormat::~LPMDFormat() { delete linecounter; }

void LPMDFormat::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la lectura/escritura de archivos en formato  \n";
 std::cout << " lpmd, este es un formato con posiciones escaladas y propio de lpmd.           \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso lpmd.                                       \n";
 std::cout << "      file          : Especifica el archivo que posee el formato lpmd.         \n";
 std::cout << "      level         : Se especifica el nivel del formato de lpmd, estos son    \n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=lpmd file=inputfile.lpmd level=0                                 \n";
 std::cout << " output module=lpmd file=outputfile.lpmd level=1 each=5                      \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato lpmd, en el    \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
}

void LPMDFormat::ReadHeader(std::istream & is) const
{
 std::string tmp;
 getline(is, tmp);
 (*linecounter) = 1;
 if (tmp.substr(0, 5) != "LPMD ") throw PluginError("lpmd", "File "+readfile+" doesn't seem to be in LPMD format (wrong header)");
}

// 
// Reads a configuration from a LPMD file 
//
bool LPMDFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 std::string tmp;
 getline(is, tmp);                                     // Numero de atomos
 (*linecounter)++;
 std::vector<std::string> words = SplitTextLine(tmp); 
 if (words.size() == 0) return false;
 long int natoms = atoi(words[0].c_str());
 getline(is, tmp);                                     // Vectores de la celda
 (*linecounter)++;
 words = SplitTextLine(tmp); 
 if(words.size()==9)
 {
  if (GetString("replacecell") == "true")
  {
   sc.SetVector(0, Vector(atof(words[0].c_str()), atof(words[1].c_str()), atof(words[2].c_str())));
   sc.SetVector(1, Vector(atof(words[3].c_str()), atof(words[4].c_str()), atof(words[5].c_str())));
   sc.SetVector(2, Vector(atof(words[6].c_str()), atof(words[7].c_str()), atof(words[8].c_str())));
  }
 }
 else if(words.size()==6)
 {
  if (GetString("replacecell") == "true")
  {
   sc.ReSet(atof(words[0].c_str()),atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str())*M_PI/180,atof(words[4].c_str())*M_PI/180,atof(words[5].c_str())*M_PI/180);
  }
 }
 else throw PluginError("lpmd", "Error ocurred when reading the base vectors, file \""+readfile+"\", line "+ToString<int>(*linecounter));
 long int atomcount = 0;
 for (long int i=0;i<natoms;++i)
 {
  getline(is, tmp);
  (*linecounter)++;
  words = SplitTextLine(tmp); 
  if (words.size() == 0) { }
  else if (words.size() == 4 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   sc.Create(new Atom(N));
   sc.SetFracPosition(atomcount++, pos);
   //Falta Asignar aca la propiedad del atomo.
  }
  else if(words.size() == 7 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   sc.Create(new Atom(N));
   sc.SetFracPosition(atomcount, pos);
   sc.SetVelocity(atomcount++, vel);
   //Falta asignar la propiedad del atomo.
  }
  else if(words.size() == 10 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Vector ace(atof(words[7].c_str()),atof(words[8].c_str()),atof(words[9].c_str()));
   sc.Create(new Atom(N));
   sc.SetFracPosition(atomcount, pos);
   sc.SetVelocity(atomcount, vel);
   sc.SetAcceleration(atomcount++, ace);
   //Falta asignar la propiedad del atomo.
  }
  else throw PluginError("lpmd", "An unidentified line was found in the file \""+readfile+"\", line "+ToString<int>(*linecounter));
 }
 return true;
}

void LPMDFormat::WriteHeader(std::ostream & os, std::vector<SimulationCell> * cells) const
{
 os << "LPMD 1.0" << std::endl;
}

void LPMDFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 out << sc.size() << std::endl;
 out << sc.GetVector(0) << " " << sc.GetVector(1) << " " << sc.GetVector(2) << std::endl;
 if(level == 0)
 {
  for (unsigned long int i=0;i<sc.size();i++) 
     out << sc[i].Symb() << " " << sc.FracPosition(i) << std::endl; //<< " " << (sc.GetAtom(i)).Type() << std::endl;
 }
 else if(level == 1)
 {
  for (unsigned long int i=0;i<sc.size();i++) 
     out << sc[i].Symb() << " " << sc.FracPosition(i) << " " << sc[i].Velocity() << std::endl; //<<" "<<(sc.GetAtom(i)).Type()<< std::endl;
 }
 else if(level == 2)
 {
  for (unsigned long int i=0;i<sc.size();i++) 
     out << sc[i].Symb() << " " << sc.FracPosition(i) << " " << sc[i].Velocity() << " " << sc[i].Acceleration() << std::endl;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new LPMDFormat(args); }
void destroy(Module * m) { delete m; }

