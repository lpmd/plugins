//
//
//

#include "lpmd.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>
#include <lpmd/atom.h>

using namespace lpmd;

LPMDFormat::LPMDFormat(std::string args): Module("lpmd")
{
 AssignParameter("version", "2.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 linecounter = new long int;
 DefineKeyword("file");
 DefineKeyword("each", "1");
 DefineKeyword("level", "0");
 DefineKeyword("extra");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 rcell = GetBool("replacecell");
 extra = SplitTextLine(GetString("extra"),',');
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
 std::cout << "      extra         : Informacion extra en el fichero, valores soportados son  \n";
 std::cout << "                      colors,type.                                             \n";
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
 if (tmp.substr(0, 5) != "LPMD ") throw PluginError("lpmd", "File "+readfile+" doesn't seem to be in LPMD X.X format (wrong header)");
 if (tmp.substr(6, 3) =="1.0")
 {
  //assume 1.0 format
  std::string info;
  getline(is, info);
  getline(is, info);
  getline(is, info);
  std::vector<std::string> words = SplitTextLine(info);
  is.unget();is.unget();
  is.unget();is.unget();
  is.unget();is.unget();
  if (words.size()==4)
  {
   hdr.push_back(std::string("SYM"));
   hdr.push_back(std::string("X"));hdr.push_back(std::string("Y"));hdr.push_back(std::string("Z"));
  }
  else if (words.size()==7)
  {
   hdr.push_back(std::string("SYM"));
   hdr.push_back(std::string("X"));hdr.push_back(std::string("Y"));hdr.push_back(std::string("Z"));
   hdr.push_back(std::string("VX"));hdr.push_back(std::string("VY"));hdr.push_back(std::string("VZ"));
  }
  else if (words.size()==10)
  {
   hdr.push_back(std::string("SYM"));
   hdr.push_back(std::string("X"));hdr.push_back(std::string("Y"));hdr.push_back(std::string("Z"));
   hdr.push_back(std::string("VX"));hdr.push_back(std::string("VY"));hdr.push_back(std::string("VZ"));
   hdr.push_back(std::string("FX"));hdr.push_back(std::string("FY"));hdr.push_back(std::string("FZ"));
  }
  else
  {
   throw PluginError("lpmd", "File "+readfile+" not have a apropiate 1.0 version");
  }
 }
 else if (tmp.substr(6, 3)=="2.0")
 {
  getline(is, tmp);
  std::string info = tmp ;
  (*linecounter)++;
  if (tmp.substr(0, 4) != "HDR ") throw PluginError("lpmd", "File"+readfile+" doesn't seem to be in LPMD 2.0 fromat (wrong HDR)");
  std::vector <std::string> words = SplitTextLine(info);
  for (unsigned long int i=0;i<words.size() ; ++i)
  {
   hdr.push_back(std::string(words[i]));
  }
 }
 else 
 {
  throw PluginError("lpmd", "The level of the file "+readfile+" are not supporten in this version of lpmd plugin.");
 }
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
  else if (words.size()+1 != hdr.size())
  {
   throw PluginError("lpmd", "Error ocurred, the header and atom information not match!");
  }
  else if (words.size()+1 == hdr.size())
  {
   int N=0;
   double X=0.0e0,Y=0.0e0,Z=0.0e0;
   double VX=0.0e0,VY=0.0e0,VZ=0.0e0;
   double AX=0.0e0,AY=0.0e0,AZ=0.0e0;
   lpmd::Vector color(0,0,0);
   for (unsigned long int j=1 ; j<hdr.size() ; ++j) // note start in one because 0 is HDR.
   {
    if (hdr[j] == "SYM") {N=ElemNum(words[j]); color = GetSpcColor(N);}
    if (hdr[j] == "X") X=atof(words[j].c_str());
    if (hdr[j] == "Y") Y=atof(words[j].c_str());
    if (hdr[j] == "Z") Z=atof(words[j].c_str());
    if (hdr[j] == "VX") VX=atof(words[j].c_str());
    if (hdr[j] == "VY") VY=atof(words[j].c_str());
    if (hdr[j] == "VZ") VZ=atof(words[j].c_str());
    if (hdr[j] == "AX") AX=atof(words[j].c_str());
    if (hdr[j] == "AY") AY=atof(words[j].c_str());
    if (hdr[j] == "AZ") AZ=atof(words[j].c_str());
    if (hdr[j] == "RGB") color = Vector(words[j]);
   }
   Vector pos(X,Y,Z);
   Vector vel(VX,VY,VZ);
   Vector ace(AX,AY,AZ);
   sc.Create(new Atom(N));
   sc.SetFracPosition(atomcount, pos);
   sc.SetVelocity(atomcount, vel);
   sc.SetColor(atomcount, color);
   sc.SetAcceleration(atomcount++, ace);
   //NOTE : falta asignar la propiedad del atomo.
  }
  else throw PluginError("lpmd", "An unidentified line was found in the file \""+readfile+"\", line "+ToString<int>(*linecounter));
 }
 return true;
}

void LPMDFormat::WriteHeader(std::ostream & os, std::vector<SimulationCell> * cells) const
{
 os << "LPMD 2.0" << std::endl;
 os << "HDR ";
 if(hdr.size()<2)
 {
  //hdr not set, using the plugin information.
  hdr.clear();
  hdr.push_back("SYM");
  hdr.push_back("X");hdr.push_back("Y");hdr.push_back("Z");
  if (level>=1)
  {
   hdr.push_back("VX");hdr.push_back("VY");hdr.push_back("VZ");
  }
  if (level>=2)
  {
   hdr.push_back("AX");hdr.push_back("AY");hdr.push_back("AZ");
  }
  if (extra.size()>0)
  {
   for(unsigned int i=0;i<extra.size();++i)
   {
    hdr.push_back(extra[i].c_str());
   }
  }
 }
 for (unsigned int i=0 ; i < hdr.size() ; ++i)
 {
  os << hdr[i] << " ";
 }
 os << '\n';
}

void LPMDFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 out << sc.size() << std::endl;
 out << sc.GetVector(0) << " " << sc.GetVector(1) << " " << sc.GetVector(2) << std::endl;
 for (unsigned long int i=0;i<sc.size();i++)
 {
  if (level>=0)
  {
   out << sc[i].Symb() << " " << sc.FracPosition(i) ;
  }
  if (level>=1)
  {
   out << " "<< sc[i].Velocity();
  }
  if (level>=2)
  {
   out << " "<< sc[i].Acceleration();
  }
  if (extra.size()>=1)
  {
   for (unsigned long j=0 ; j < extra.size() ; ++j)
   {
    if(extra[j] == "color") {lpmd::Vector tmp=sc[i].Color(); out << "          "<< tmp.Write(); }
    if(extra[j] == "type") { out << "          " << "ATOMTYPE"; }
   }
  }
  out << '\n';
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new LPMDFormat(args); }
void destroy(Module * m) { delete m; }

