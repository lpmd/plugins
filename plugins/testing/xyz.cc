	//
//
//

#include "xyz.h"

#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <lpmd/manipulations.h>

using namespace lpmd;

XYZFormat::XYZFormat(std::string args): Plugin("xyz", "2.0")
{
 ParamList & params = (*this);
 linecounter = new long int;
 //
 DefineKeyword("file");
 DefineKeyword("level", "0");
 DefineKeyword("each", "1");
 DefineKeyword("coords", "positive");
 DefineKeyword("inside", "false");
 DefineKeyword("external", "consider");
 AssignParameter("external", "consider");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = params["file"];
 interval = int(params["each"]);
 level = int(params["level"]);
 coords = params["coords"];
 inside = params["inside"];
 external = params["external"];
 rcell = params["replacecell"];
}

XYZFormat::~XYZFormat() { delete linecounter; }

void XYZFormat::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la lectura/escritura de archivos en formato  \n";
 std::cout << " xyz, este es un formato con posiciones absolutas.                             \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso xyz.                                        \n";
 std::cout << "      file          : Especifica el archivo que posee el formato lpmd.         \n";
 std::cout << "      level         : Se especifica el nivel del formato de xyz, estos son     \n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << "      coords        : Especifica si la celda esta o no centrada en el          \n";
 std::cout << "                      origen (centered/positive=default).                      \n";
 std::cout << "      inside        : Especifica si se deben reacomodar los atomos que         \n";
 std::cout << "                      se encuentran fuera de la celda (true/false=default).    \n";
 std::cout << "      external      : Especifica si se deben ignorar o no los atomos que       \n";
 std::cout << "                      se encuentran fuera de la celda(ignore/consider=default).\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=xyz file=inputfile.xyz level=0                                   \n";
 std::cout << " output module=xyz file=outputfile.xyz level=1 each=5                        \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato xyz, en el     \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
}

void XYZFormat::ReadHeader(std::istream & is) const
{
 assert(&is != 0); //icc 869
 (*linecounter) = 0;
 // El formato XYZ no tiene ningun header especial
}

// 
// Lee una configuracion desde un archivo XYZ 
//
bool XYZFormat::ReadCell(std::istream & is, Configuration & sc) const
{
 if ((*this)["replacecell"] == "true") throw PluginError("xyz", "This format does not contain any cell vectors.");
 std::string tmp;

 // Tomamos las cell y los atomos.
 BasicParticleSet & part = sc.Atoms();
 BasicCell & cell = sc.Cell();
 assert(part.Size() == 0);
 //
 getline(is, tmp);              // This reads the "number of atoms" line
 (*linecounter)++;
 char * buffer = new char[100];
 char * bfold = buffer;
 buffer[0] = 'X';
 char ** endptr = &buffer;
 long natoms = strtol(tmp.c_str(), endptr, 10);
 if ((buffer[0] == '\0') && (natoms > 0)) { }
 else 
 {
  delete [] bfold;
  return false;
 }
 delete [] bfold; // hay que hacer delete sobre el antiguo valor, no sobre el nuevo ya que strtol lo cambia

 //
 getline(is, tmp);             // This reads the "title" line and ignores it
 (*linecounter)++;

 for (long i=0;i<natoms;++i)
 { 
  getline(is, tmp);
  (*linecounter)++;
  Array<std::string> words = StringSplit(tmp, ' '); 
  if (words.Size() == 4)
  {
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   part.Append(Atom(words[0], pos));
   sc.SetTag(sc, Tag("level"), 0);
  }
  else if (words.Size() == 7)
  {
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   part.Append(Atom(words[0], pos, vel));
   sc.SetTag(sc, Tag("level"), 1);
  }
  else if (words.Size() == 10)
  {
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Vector ace(atof(words[7].c_str()),atof(words[8].c_str()),atof(words[9].c_str()));
   part.Append(Atom(words[0], pos, vel, ace));
   sc.SetTag(sc, Tag("level"), 2);
  }
  else throw PluginError("xyz", "An unidentified line was found in the file \""+readfile+"\", line "+ToString<int>(*linecounter));
 }
 if (coords == "centered") UnCenter(part,cell);
 if (external != "ignore")
 {
  if (inside == "true")
  {
   //Reubica los atomos dentro de la celda.
   for(int i=0;i<natoms;i++)
   {
    part[i].Position() = cell.FittedInside(part[i].Position());
   }
  }
  //Chequea que todos los atomos que se han leido esten dentro de la "celda".
  else
  {
   for (int i=0;i<natoms;i++)
   {
    if (!cell.IsInside(part[i].Position()))
    {
     throw PluginError("xyz", "The atom["+ToString<int>(i)+"] was found outside the cell.");
    }
   }
  }
 }
 return true;
}

void XYZFormat::WriteHeader(std::ostream & os, SimulationHistory * sh) const
{
 assert(&os != 0); //icc 869
 assert(sh >= (void *)NULL); //icc 869
 //assert(&sh != 0); //icc 869
 // El formato XYZ no tiene ningun header especial
}

void XYZFormat::WriteCell(std::ostream & out, Configuration & sc) const
{
 // Tomamos las cell y los atomos.
 BasicParticleSet & part = sc.Atoms();
 BasicCell & cell = sc.Cell();
 if (inside == "true")
 {
  //Reubica los atomos dentro de la celda.
  for (long int i=0;i<part.Size();i++) part[i].Position() = cell.FittedInside(part[i].Position());
 }
 out << part.Size() << std::endl;
 out << '\n';
 if (level == 0)
 {
  for(long int i=0;i<part.Size();i++) out << part[i].Symbol() << " " << part[i].Position() << std::endl;
 }
 else if (level == 1)
 {
  for (long int i=0;i<part.Size();i++) out << part[i].Symbol() << " " << part[i].Position() << " " << part[i].Velocity() << std::endl;
 }
 else if (level == 2)
 {
  for (long int i=0;i<part.Size();i++) out << part[i].Symbol() << " " << part[i].Position() << " " << part[i].Velocity() << " " << part[i].Acceleration() << std::endl;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new XYZFormat(args); }
void destroy(Plugin * m) { delete m; }

