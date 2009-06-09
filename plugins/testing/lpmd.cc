//
//
//

#include "lpmd.h"

#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <lpmd/atom.h>

using namespace lpmd;

LPMDFormat::LPMDFormat(std::string args): Plugin("lpmd", "2.0")
{
 ParamList & params = (*this);
 //
 linecounter = new long int;
 DefineKeyword("file");
 DefineKeyword("each", "1");
 DefineKeyword("level", "0");
 DefineKeyword("extra");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = (*this)["file"];
 interval = int(params["each"]);
 level = int(params["level"]);
 rcell = params["replacecell"];
 extra = StringSplit((*this)["extra"],',');
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
 std::cout << "                      RGB,C,TYPE.                                              \n";
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
 if (tmp.substr(5, 3) =="1.0")
 {
  //assume 1.0 format
  int where = is.tellg();
  std::string info;
  getline(is, info);
  getline(is, info);
  getline(is, info);
  Array<std::string> words = StringSplit(info, ' ');
  is.seekg(where);
  if (words.Size()==4)
  {
   hdr.Append(std::string("HDR"));
   hdr.Append(std::string("SYM"));
   hdr.Append(std::string("X"));hdr.Append(std::string("Y"));hdr.Append(std::string("Z"));
  }
  else if (words.Size()==7)
  {
   hdr.Append(std::string("HDR"));
   hdr.Append(std::string("SYM"));
   hdr.Append(std::string("X"));hdr.Append(std::string("Y"));hdr.Append(std::string("Z"));
   hdr.Append(std::string("VX"));hdr.Append(std::string("VY"));hdr.Append(std::string("VZ"));
  }
  else if (words.Size()==10)
  {
   hdr.Append(std::string("HDR"));
   hdr.Append(std::string("SYM"));
   hdr.Append(std::string("X"));hdr.Append(std::string("Y"));hdr.Append(std::string("Z"));
   hdr.Append(std::string("VX"));hdr.Append(std::string("VY"));hdr.Append(std::string("VZ"));
   hdr.Append(std::string("FX"));hdr.Append(std::string("FY"));hdr.Append(std::string("FZ"));
  }
  else
  {
   throw PluginError("lpmd", "File "+readfile+" not have a apropiate 1.0 version");
  }
 }
 else if (tmp.substr(5, 3)=="2.0")
 {
  getline(is, tmp);
  std::string info = tmp ;
  (*linecounter)++;
  if (tmp.substr(0, 4) != "HDR ") throw PluginError("lpmd", "File"+readfile+" doesn't seem to be in LPMD 2.0 fromat (wrong HDR)");
  Array<std::string> words = StringSplit(info,' ');
  for (long int i=0;i<words.Size() ; ++i)
  {
   hdr.Append(std::string(words[i]));
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
bool LPMDFormat::ReadCell(std::istream & is, Configuration & con) const
{
 BasicCell & cell = con.Cell();
 BasicParticleSet & part = con.Atoms();
 assert(part.Size() == 0);
 std::string tmp;
 getline(is, tmp);                                     // Numero de atomos
 (*linecounter)++;
 Array<std::string> words = StringSplit(tmp, ' '); 
 if (words.Size() == 0) return false;
 long int natoms = atoi(words[0].c_str());
 getline(is, tmp);                                     // Vectores de la celda
 (*linecounter)++;
 words = StringSplit(tmp, ' '); 
 if(words.Size()==9)
 {
  if ((*this)["replacecell"] == "true")
  {
   cell[0] = Vector(atof(words[0].c_str()), atof(words[1].c_str()), atof(words[2].c_str()));
   cell[1] = Vector(atof(words[3].c_str()), atof(words[4].c_str()), atof(words[5].c_str()));
   cell[2] = Vector(atof(words[6].c_str()), atof(words[7].c_str()), atof(words[8].c_str()));
  }
 }
 else if(words.Size()==6)
 {
  if ((*this)["replacecell"] == "true")
  {
   Cell tmp(atof(words[0].c_str()),atof(words[1].c_str()),atof(words[2].c_str()),atof(words[3].c_str())*M_PI/180,atof(words[4].c_str())*M_PI/180,atof(words[5].c_str())*M_PI/180);
   for (int q=0;q<3;++q) cell[q] = tmp[q];
  }
 }
 else throw PluginError("lpmd", "Error ocurred when reading the base vectors, file \""+readfile+"\", line "+ToString<int>(*linecounter));
 long int atomcount = 0;
 for (long int i=0;i<natoms;++i)
 {
  getline(is, tmp);
  (*linecounter)++;
  words = StringSplit(tmp, ' ');
  if (words.Size() == 0) 
  {
   throw PluginError("lpmd", "Error ocurred, the atom file not have elements!"); 
  }
  else if (words.Size() >=1)
  {
   std::string sym;
   double X=0.0e0,Y=0.0e0,Z=0.0e0;
   double VX=0.0e0,VY=0.0e0,VZ=0.0e0;
   double AX=0.0e0,AY=0.0e0,AZ=0.0e0;
   lpmd::Vector color(0,0,0);
   double colors=-1.0e0;
   for (long int k=1 ; k < hdr.Size() ; ++k)
   {
#warning Color seteado a cero mientras tanto, para no usar GetSpcColor
    if (hdr[k] == "SYM") {sym=words[k-1]; color = Vector(0,0,0); } //GetSpcColor(N);}
    if (hdr[k] == "X") X=atof(words[k-1].c_str());
    if (hdr[k] == "Y") Y=atof(words[k-1].c_str());
    if (hdr[k] == "Z") Z=atof(words[k-1].c_str());
    if (hdr[k] == "VX") VX=atof(words[k-1].c_str());
    if (hdr[k] == "VY") VY=atof(words[k-1].c_str());
    if (hdr[k] == "VZ") VZ=atof(words[k-1].c_str());
    if (hdr[k] == "AX") AX=atof(words[k-1].c_str());
    if (hdr[k] == "AY") AY=atof(words[k-1].c_str());
    if (hdr[k] == "AZ") AZ=atof(words[k-1].c_str());
    if (hdr[k] == "RGB") color = Vector(words[k-1].c_str());
    if (hdr[k] == "C") colors=atof(words[k-1].c_str());
   }
   Vector pos = cell.Cartesian(Vector(X,Y,Z));
   Vector vel(VX,VY,VZ);
   Vector ace(AX,AY,AZ);
   part.Append(Atom(sym,pos,vel,ace));
#warning Hay que asignar el color de alguna forma!!!
   atomcount++;
   //part.SetFracPosition(atomcount, pos);
   //part.SetVelocity(atomcount, vel);
   //part.SetColor(atomcount, color);
   //if(colors>=0 && colors <=1) sc[i].SetColor(colors);
   //sc.SetAcceleration(atomcount++, ace);
#warning Falta asignar la propiedad del atomo?
  }
  else throw PluginError("lpmd", "An unidentified line was found in the file \""+readfile+"\", line "+ToString<int>(*linecounter));
 }
 return true;
}

void LPMDFormat::WriteHeader(std::ostream & os, SimulationHistory * sh) const
{
 os << "LPMD 2.0" << std::endl;
 os << "HDR ";
 if(hdr.Size()<2)
 {
  //hdr not set, using the plugin information.
  hdr.Clear();
  hdr.Append("SYM");
  hdr.Append("X");hdr.Append("Y");hdr.Append("Z");
  if (level>=1)
  {
   hdr.Append("VX");hdr.Append("VY");hdr.Append("VZ");
  }
  if (level>=2)
  {
   hdr.Append("AX");hdr.Append("AY");hdr.Append("AZ");
  }
  if (extra.Size()>0)
  {
   for(long int i=0;i<extra.Size();++i)
   {
    hdr.Append(extra[i].c_str());
   }
  }
 }
 for (long int i=0 ; i < hdr.Size() ; ++i)
 {
  os << hdr[i] << " ";
 }
 os << '\n';
}

void LPMDFormat::WriteCell(std::ostream & out, Configuration & con) const
{
 BasicParticleSet & part = con.Atoms();
 BasicCell & cell = con.Cell();

 level = int(Parameter(con.GetTag(con,"level")));
 out << part.Size() << std::endl;
 out << cell[0] << " " << cell[1] << " " << cell[2] << std::endl;
 for (long int i=0;i<part.Size();i++)
 {
  if (level>=0)
  {
   out << part[i].Symbol() << " " << cell.Fractional(part[i].Position()) ;
  }
  if (level>=1)
  {
   out << " "<< part[i].Velocity();
  }
  if (level>=2)
  {
   out << " "<< part[i].Acceleration();
  }
  if (extra.Size()>=1)
  {
   for (long int j=0 ; j < extra.Size() ; ++j)
   {
#warning Comentadas Opciones COLOR para lpmd2
    //if(extra[j] == "RGB") {lpmd::Vector tmp=part[i].Color(); FormattedWrite(out,tmp); }
    //if(extra[j] == "C") {double color=part[i].ColorS(); out << "          " << color;} 
    if(extra[j] == "TYPE") { out << "          " << "ATOMTYPE"; }
   }
  }
  out << '\n';
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LPMDFormat(args); }
void destroy(Plugin * m) { delete m; }

