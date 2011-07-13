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
 DefineKeyword("zerocm", "false");
 DefineKeyword("external", "consider");
 DefineKeyword("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = params["file"];
 interval = int(params["each"]);
 level = int(params["level"]);
 coords = params["coords"];
 inside = params["inside"];
 external = params["external"];
 rcell = params["replacecell"];
 zerocm = params["zerocm"];
}

XYZFormat::~XYZFormat() { delete linecounter; }

void XYZFormat::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = xyz                                                      \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugins is used for read/write files in xyz format, this format have\n";
 std::cout << " the atomic positions in [A].                                                  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      file          : Specify the filename of the xyz input.                   \n";
 std::cout << "      level         : Set the level of the xyz file, the levels availables are \n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-acce.                      \n";
 std::cout << "      coords        : Indicate if the input cell is centered in the origin or  \n";
 std::cout << "                      not, the possible values are centered/positive=default.  \n";
 std::cout << "      inside        : Specify if the atoms outside of the cell boundary have to\n";
 std::cout << "                      bee rellocated inside (periodicity) the available values \n";
 std::cout << "                      are true/false=default.                                  \n";
 std::cout << "      external      : Indicate if you have to ignore or not the atoms that are \n";
 std::cout << "                      located outside of the simulationc cell, the availables  \n";
 std::cout << "                      values are ignore/consider=default.                      \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading a inputfile using this plugin.                                       \n";
 std::cout << " input module=xyz file=inputfile.xyz level=0                                   \n";
 std::cout << " #Writing using this plugin.                                                   \n";
 std::cout << " output module=xyz file=outputfile.xyz level=1 each=5                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 if (zerocm == "true")
 {
  Vector cm(0.0, 0.0, 0.0);
  const Vector geocenter = cell.Cartesian(Vector(0.5, 0.5, 0.5));
  double totalmass = 0.0;
  for (long i=0;i<natoms;++i) 
  {
   cm += (part[i].Mass()*part[i].Position());
   totalmass += part[i].Mass();
  }
  for (long i=0;i<natoms;++i) part[i].Position() += (geocenter-(cm/totalmass));
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


bool XYZFormat::SkipCell(std::istream & is) const
{
 std::string tmp;
 getline(is, tmp);              // This reads the "number of atoms" line
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
 getline(is, tmp);
 for (long i=0;i<natoms;++i) getline(is, tmp);
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

