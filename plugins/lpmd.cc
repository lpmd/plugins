//
//
//

#include <lpmd/util.h>

#include "lpmd.h"

using namespace lpmd;

LPMDFormat::LPMDFormat(std::string args): Module("lpmd")
{
 interval = 1;
 level = 0;
 ProcessArguments(args);
}

LPMDFormat::~LPMDFormat() { }

void LPMDFormat::SetParameter(std::string name)
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
}

void LPMDFormat::Show() const
{
 Module::Show();
 std::cout << "   file     = " << writefile << '\n';
 std::cout << "   each     = " << interval << '\n';
 std::cout << "   level    = " << level << '\n';
}

void LPMDFormat::ShowHelp() const
{
}

std::string LPMDFormat::Keywords() const 
{
 return "file each level";
}

void LPMDFormat::ReadHeader(std::istream & is) const
{
 std::string tmp;
 getline(is, tmp);
 // FIXME: Comparar aqui la linea de encabezado para ver si comienza con "LPMD" y checkear numero de version
}

// 
// Reads a configuration from a LPMD file 
//
void LPMDFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 std::string tmp;
 getline(is, tmp);                                     // Numero de atomos
 std::vector<std::string> words = SplitTextLine(tmp); 
 sc.Initialize(atoi(words[0].c_str()));                // Clear and Reserve space for the atoms
 getline(is, tmp);                                     // Vectores de la celda
 words = SplitTextLine(tmp); 
 if(words.size()==9)
 {
  sc.SetVector(0, Vector(atof(words[0].c_str()), atof(words[1].c_str()), atof(words[2].c_str())));
  sc.SetVector(1, Vector(atof(words[3].c_str()), atof(words[4].c_str()), atof(words[5].c_str())));
  sc.SetVector(2, Vector(atof(words[6].c_str()), atof(words[7].c_str()), atof(words[8].c_str())));
 }
 else if(words.size()==6)
 {
  sc.ReSet(atof(words[0].c_str()),atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str())*M_PI/180,atof(words[4].c_str())*M_PI/180,atof(words[5].c_str())*M_PI/180);
 }
 else {std::cerr << "Error ocurred during LPMD format lecture in Base-Vector line"<<'\n';}
 long int atomcount = 0;
 while(getline(is, tmp))
 {
  words = SplitTextLine(tmp); 
  if (words.size() == 0) { }
  else if (words.size() == 4 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   sc.AppendAtom(Atom(N));
   sc.SetFracPosition(atomcount++, pos);
   //Falta Asignar aca la propiedad del atomo.
  }
  else if(words.size() == 7 && words[0] != "#")
  {
   int N=ElemNum(words[0]);
   Vector pos(atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str()));
   Vector vel(atof(words[4].c_str()),atof(words[5].c_str()),atof(words[6].c_str()));
   sc.AppendAtom(Atom(N));
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
   sc.AppendAtom(Atom(N));
   sc.SetFracPosition(atomcount, pos);
   sc.SetVelocity(atomcount, vel);
   sc.SetAcceleration(atomcount++, ace);
   //Falta asignar la propiedad del atomo.
  }
  else EndWithError("A unidentified line was found in the input file.");
 }
}

void LPMDFormat::WriteHeader(std::ostream & os) const
{
 os << "LPMD 1.0" << std::endl;
}

void LPMDFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 out << sc.Size() << std::endl;
 out << sc.GetVector(0) << " " << sc.GetVector(1) << " " << sc.GetVector(2) << std::endl;
 if(level == 0)
 {
  for(int i=0;i<sc.Size();i++) out <<(sc.GetAtom(i)).Symb()<<" "<< (sc.FracPosition(i)) << std::endl; //<< " " << (sc.GetAtom(i)).Type() << std::endl;
 }
 else if(level == 1)
 {
  for (int i=0;i<sc.Size();i++) out <<(sc.GetAtom(i)).Symb()<<" "<< (sc.FracPosition(i)) << " " << (sc.GetAtom(i)).Velocity() << std::endl; //<<" "<<(sc.GetAtom(i)).Type()<< std::endl;
 }
 else if(level == 2)
 {
  for(int i=0;i<sc.Size();i++) out <<(sc.GetAtom(i)).Symb()<< (sc.FracPosition(i)) << " " << (sc.GetAtom(i)).Velocity() << " " << (sc.GetAtom(i)).Acceleration() << std::endl;
 }
}


// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new LPMDFormat(args); }
void destroy(Module * m) { delete m; }


