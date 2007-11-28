//
//
//

#include <lpmd/util.h>

#include "xyz.h"

using namespace lpmd;

XYZFormat::XYZFormat(std::string args): Module("xyz")
{
 interval = 1;
 level = 0;
 coords = "positive";
 ProcessArguments(args);
}

XYZFormat::~XYZFormat() { }

void XYZFormat::SetParameter(std::string name)
{
 if (name == "file")
 {
  AssignParameter("file", GetNextWord());
  readfile = GetString("file");
  writefile = readfile;
 }
 else if (name == "each")
 {
  AssignParameter("each", GetNextWord());
  interval = GetInteger("each");
 }
 else if (name == "level")
 {
  AssignParameter("level", GetNextWord());
  level = GetInteger("level");
 }
 else if (name == "coords")
 {
  AssignParameter("coords", GetNextWord());
  coords = GetString("coords");
 }
}

void XYZFormat::Show() const
{
 Module::Show();
 std::cout << "   file     = " << writefile << '\n';
 if (Defined("interval")) std::cout << "   interval = " << interval << '\n';
 if (Defined("level")) std::cout << "   level    = " << level << '\n';
}

void XYZFormat::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = XYZ " << '\n';
 std::cout << " Module Version     = 1.0 " << '\n';
 std::cout << " Support API lpmd   = 1.0 " << '\n';
 std::cout << " Problems Report to = gnm@gnm.cl " << '\n';
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Information "<< '\n' << '\n';
 std::cout << "      El modulo ayuda a lpmd leer y escribir ficheros xyz. \n";
 std::cout << " para ello lpmd necesita saber el NIVEL del formato para \n";
 std::cout << " poder leer o escribir correctamente. " << '\n';
 std::cout << "      Los niveles de XYZ son 0,1 y 2 donde 0 muestra posiciones \n";
 std::cout << " 1 muestra posiciones y velocidades, 2 incluye aceleraciones\n\n";
 std::cout << " Use to read    = input module=xyz file=file-input.xyz level=0/1/2 \n" ;
 std::cout << " Use to write   = output module=xyz file=file-output.xyz each=steps level=0/1/2\n" ;
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example " << '\n' << '\n';
 std::cout << " Escribiendo    = output module=xyz file=fileoutput.xyz each=10 level=0\n" << '\n';
 std::cout << " En este caso el fichero de salida es de la forma xyz y se graba\n";
 std::cout << " cada 10 steps con un nivel 0" <<'\n';
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string XYZFormat::Keywords() const
{
 return "file each level coords";
}

void XYZFormat::ReadHeader(std::istream & is) const
{
 // El formato XYZ no tiene ningun header especial
}

// 
// Lee una configuracion desde un archivo XYZ 
//
void XYZFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 std::string tmp;
 while(getline(is, tmp))
 {
  std::vector<std::string> words = SplitTextLine(tmp); 
  if (words.size() == 1 && words[0]!="#") sc.Initialize(atoi(words[0].c_str()));      // Clear and Reserve space for the atoms
  else if (words.size() == 0) { }
  else if (words.size() == 4 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Atom a(N,pos);
   sc.AppendAtom(a);
  }
  else if(words.size() == 7 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Atom a(N,pos,vel);
   sc.AppendAtom(a);
  }
  else if(words.size() == 10 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   Vector ace(atof(words[7].c_str()),atof(words[8].c_str()),atof(words[9].c_str()));
   Atom a(N,pos,vel);
   sc.AppendAtom(a);
  }
  else EndWithError("A unidentified line was found in the input file.");
 }
 if (coords == "centered") sc.UnCenter();
}

void XYZFormat::WriteHeader(std::ostream & os) const
{
 // El formato XYZ no tiene ningun header especial
}

void XYZFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
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
