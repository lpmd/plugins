	//
//
//

#include "xyz.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>

using namespace lpmd;

XYZFormat::XYZFormat(std::string args): Module("xyz")
{
 linecounter = new long int;
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
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
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 coords = GetString("coords");
 inside = GetString("inside");
 external = GetString("external");
 rcell = GetBool("replacecell");
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
 (*linecounter) = 0;
 // El formato XYZ no tiene ningun header especial
}

// 
// Lee una configuracion desde un archivo XYZ 
//
bool XYZFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 if (GetString("replacecell") == "true") throw PluginError("xyz", "This format does not contain any cell vectors.");
 std::string tmp;

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

 // 
 for (long i=0;i<natoms;++i)
 { 
  getline(is, tmp);
  (*linecounter)++;
  std::vector<std::string> words = SplitTextLine(tmp); 
  if (words.size() == 4)
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   sc.Create(new Atom(N, pos));
  }
  else if (words.size() == 7)
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   sc.Create(new Atom(N, pos, vel));
  }
  else if (words.size() == 10)
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Vector ace(atof(words[7].c_str()),atof(words[8].c_str()),atof(words[9].c_str()));
   sc.Create(new Atom(N, pos, vel, ace));
  }
  else throw PluginError("xyz", "An unidentified line was found in the file \""+readfile+"\", line "+ToString<int>(*linecounter));
 }
 if (coords == "centered") sc.UnCenter();
 if (external != "ignore")
 {
  if (inside == "true")
  {
   //Reubica los atomos dentro de la celda.
   for(int i=0;i<natoms;i++)
   {
    Vector atompos = sc[i].Position();
    sc.SetPosition(i,atompos);
   }
  }
  //Chequea que todos los atomos que se han leido esten dentro de la "celda".
  else
  {
   for (int i=0;i<natoms;i++)
   {
    Vector pos = sc[i].Position();
    sc.ConvertToInternal(pos);
    Vector opos = pos;
    for (int q=0;q<3;++q) pos[q] = pos[q]/sc.GetVector(q).Module();
    if (pos[0]<0 || pos[0]>1.0) throw PluginError("xyz", "The atom ["+ToString<int>(i)+"] was found outside the cell in the [a] Vector");
    if (pos[1]<0 || pos[1]>1.0) throw PluginError("xyz", "The atom ["+ToString<int>(i)+"] was found outside the cell in the [b] Vector");
    if (pos[2]<0 || pos[2]>1.0) throw PluginError("xyz", "The atom ["+ToString<int>(i)+"] was found outside the cell in the [c] Vector");
   }
  }
 }
 return true;
}

void XYZFormat::WriteHeader(std::ostream & os, std::vector<SimulationCell> * cells) const
{
 // El formato XYZ no tiene ningun header especial
}

void XYZFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 if (inside == "true")
 {
  //Reubica los atomos dentro de la celda.
  for (unsigned long int i=0;i<sc.size();i++) sc.SetPosition(i, sc[i].Position());
 }
 out << sc.size() << std::endl;
 out << '\n';
 if (level == 0)
 {
  for(unsigned long int i=0;i<sc.size();i++) out << sc[i].Symb() << " " << sc[i].Position() << std::endl;
 }
 else if (level == 1)
 {
  for (unsigned long int i=0;i<sc.size();i++) out << sc[i].Symb() << " " << sc[i].Position() << " " << sc[i].Velocity() << std::endl;
 }
 else if (level == 2)
 {
  for (unsigned long int i=0;i<sc.size();i++) out << sc[i].Symb() << " " << sc[i].Position() << " " << sc[i].Velocity() << " " << sc[i].Acceleration() << std::endl;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new XYZFormat(args); }
void destroy(Module * m) { delete m; }

