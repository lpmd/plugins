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
 AssignParameter("level", "0");
 AssignParameter("each", "1");
 AssignParameter("coords", "positive");
 AssignParameter("inside", "false");
 AssignParameter("external", "consider");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 coords = GetString("coords");
 inside = GetString("inside");
 external = GetString("external");
}

XYZFormat::~XYZFormat() { delete linecounter; }

void XYZFormat::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = xyz                                                      \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "                      origen (centered/uncentered=default).                    \n";
 std::cout << "      inside        : Especifica si se deben reacomodar los atomos que         \n";
 std::cout << "                      se encuentran fuera de la celda (true/false=default).    \n";
 std::cout << "      external      : Especifica si se deben ignorar o no los atomos que       \n";
 std::cout << "                      se encuentran fuera de la celda(ignore/consider=default).\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=xyz file=inputfile.xyz level=0                                   \n";
 std::cout << " output module=xyz file=outputfile.xyz level=1 each=5                        \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato xyz, en el     \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

}

std::string XYZFormat::Keywords() const
{
 return "file each level coords inside external";
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
 std::string tmp;

 //
 getline(is, tmp);              // This reads the "number of atoms" line
 (*linecounter)++;
 char * buffer = new char[100];
 char * bfold = buffer;
 buffer[0] = 'X';
 char ** endptr = &buffer;
 long natoms = strtol(tmp.c_str(), endptr, 10);
 if ((buffer[0] == '\0') && (natoms > 0)) sc.Initialize(natoms);           // Clear and Reserve space for the atoms
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
   Atom a(N,pos);
   sc.AppendAtom(a);
  }
  else if (words.size() == 7)
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Atom a(N,pos,vel);
   sc.AppendAtom(a);
  }
  else if (words.size() == 10)
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Vector ace(atof(words[7].c_str()),atof(words[8].c_str()),atof(words[9].c_str()));
   Atom a(N,pos,vel);
   sc.AppendAtom(a);
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
   double lx = (sc.GetVector(0)).Mod();
   double ly = (sc.GetVector(1)).Mod();
   double lz = (sc.GetVector(2)).Mod();
   for(int i=0;i<natoms;i++)
   {
    Vector pos=sc[i].Position();
    if (pos.GetX()<0 || pos.GetX()>lx) throw PluginError("xyz", "The atom ["+ToString<int>(i)+"] was found outside the cell in the [a] Vector");
    if (pos.GetY()<0 || pos.GetY()>ly) throw PluginError("xyz", "The atom ["+ToString<int>(i)+"] was found outside the cell in the [b] Vector");
    if (pos.GetZ()<0 || pos.GetZ()>lz) throw PluginError("xyz", "The atom ["+ToString<int>(i)+"] was found outside the cell in the [c] Vector");
   }
  }
 }
 return true;
}

void XYZFormat::WriteHeader(std::ostream & os) const
{
 // El formato XYZ no tiene ningun header especial
}

void XYZFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 if (inside == "true")
 {
  //Reubica los atomos dentro de la celda.
  for(int i=0;i<sc.Size();i++)
  {
   Vector atompos = sc[i].Position();
   sc.SetPosition(i,atompos);
  }
 }
 out << sc.Size() << std::endl;
 out << '\n';
 if(level == 0)
 {
  for(int i=0;i<sc.Size();i++) out << (sc.GetAtom(i)).Symb() << " " << (sc.GetAtom(i)).Position() << std::endl;
 }
 else if(level == 1)
 {
  for(int i=0;i<sc.Size();i++) out << (sc.GetAtom(i)).Symb() << " " << (sc.GetAtom(i)).Position() << " " << (sc.GetAtom(i)).Velocity() << std::endl;
 }
 else if(level == 2)
 {
  for(int i=0;i<sc.Size();i++) out << (sc.GetAtom(i)).Symb() << " " << (sc.GetAtom(i)).Position() << " " << (sc.GetAtom(i)).Velocity() << " " << (sc.GetAtom(i)).Acceleration() << std::endl;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new XYZFormat(args); }
void destroy(Module * m) { delete m; }

