//
//
//

#include "dlpoly.h"

#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <sstream>
#include <iomanip>
#include <lpmd/simulationhistory.h>

using namespace lpmd;

DlPolyFormat::DlPolyFormat(std::string args): Plugin("dlpoly", "2.0"), dt(0.0)
{
 lpmd::ParamList & param = (*this);
 AssignParameter("each", "1");
 AssignParameter("level", "0");
 AssignParameter("periodicity", "1");
 AssignParameter("replacecell", "false");
 AssignParameter("ftype","CONFIG");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = param["file"];
 interval = int(param["each"]);
 level = int(param["level"]);
 pbkey = int(param["periodicity"]);
 rcell = bool(param["replacecell"]);
 ftype = param["ftype"];
 if (Defined("dt")) dt = double(param["dt"]);
 stepcnt = new long int;
}

DlPolyFormat::~DlPolyFormat() { delete stepcnt; }

void DlPolyFormat::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = dlpoly                                                   \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la lectura/escritura de archivos en formato  \n";
 std::cout << " CONFIG de DLPOLY, este es un formato con posiciones absolutas.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      file          : Especifica el archivo que posee el formato dlpoly.       \n";
 std::cout << "      level         : Se especifica el nivel del formato del archivo, estos son\n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << "      periodicity   : Flag periodic boundary key del formato CONFIG            \n";
 std::cout << "      ftype         : File type, puede ser CONFIG o HISTORY.                   \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=dlpoly file=CONFIG level=0                                       \n";
 std::cout << " output module=dlpoly file=CONFIG.output level=1 periodicity=2 each=5          \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato DLPOLY, en el  \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void DlPolyFormat::ReadHeader(std::istream & is) const
{
 // El formato CONFIG de DLPOLY no tiene ningun header especial
 // El formato HISTORY de DLPOLY posee un titulo y un set de valores de la simulacion.
 if(ftype=="HISTORY")
 {
  std::string tmp;
  getline(is,tmp); // lee la linea de titulo.
  getline(is,tmp);
  std::istringstream ost(tmp);
  int a1,a2,a3,a4,a5;
  ost >> a1 >> a2 >> a3 >> a4 >> a5;
  if (a1!=level) throw PluginError("dlpoly", "Error the level-file are diferent that the specified.");
 }
}

// 
// Lee una configuracion desde un archivo CONFIG 
//
bool DlPolyFormat::ReadCell(std::istream & is, Configuration & con) const
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 con.SetTag(con, Tag("level"), level);
 int fkey, pbk;
 std::string tmp;
 Vector cv[3];
 double x, y, z;
 if(ftype=="CONFIG")
 {
  getline(is, tmp);           // lee la linea de titulo
  if (is.eof()) return false; // no hay mas configuraciones que leer
  getline(is, tmp);
  std::istringstream ost(tmp);
  ost >> fkey >> pbk;
  for (int i=0;i<3;++i)
  {
   getline(is, tmp);
   std::istringstream vst(tmp);
   vst >> x >> y >> z;
   cv[i] = Vector(x, y, z);
   if ((*this)["replacecell"] == "true") cell[i] = cv[i];
  } 
  long cc = 0;
  while (1)
  {
   cc++;
   getline(is, tmp);
   if (is.eof()) break;
   std::istringstream vst(tmp); 
   std::string symbol;
   vst >> symbol;
   Atom this_atom(ElemNum(symbol));
   for (int i=0;i<=fkey;++i)
   {
    getline(is, tmp);
    std::istringstream vst(tmp);
    vst >> x >> y >> z;
    if (i == 0) this_atom.Position() = Vector(x+0.5*cv[0].Module(), y+0.5*cv[1].Module(), z+0.5*cv[2].Module());
    if (i == 1) this_atom.Velocity() = Vector(x, y, z);
    if (i == 2) this_atom.Acceleration() = Vector(x, y, z);
   }
   atoms.Append(this_atom);
  }
 }
 else if(ftype=="HISTORY")
 {
  getline(is, tmp);
  if (is.eof()) return false; // no hay mas configuraciones que leer
  std::istringstream tst(tmp);
  std::string ttmp;
  int megatm;
  tst >> ttmp >> ttmp ;
  tst >> megatm ;
  tst >> fkey >> pbk ;
  for(int j=0;j<3;++j)
  {
   getline(is, tmp);
   std::istringstream vst(tmp);
   vst >> x >> y >> z;
   cv[j] = Vector(x,y,z);
   if((*this)["replacecell"] == "true") cell[j] = cv[j];
  }
  for(int k=0;k<megatm;++k)
  {
   getline(is,tmp);
   std::istringstream vst(tmp); 
   std::string symbol;
   vst >> symbol;
   Atom this_atom(ElemNum(symbol));
   for (int i=0;i<=fkey;++i)
   {
    getline(is, tmp);
    std::istringstream vst(tmp);
    vst >> x >> y >> z;
    if (i == 0) this_atom.Position() = Vector(x+0.5*cv[0].Module(), y+0.5*cv[1].Module(), z+0.5*cv[2].Module());
    if (i == 1) this_atom.Velocity() = Vector(x, y, z);
    if (i == 2) this_atom.Acceleration() = Vector(x, y, z);
   }
   atoms.Append(this_atom);
  }
 }
 else
 {
  throw PluginError("dlpoly", "The fkey="+ftype+" is not a valid fkey value.");
 }
 return true;
}

void DlPolyFormat::WriteHeader(std::ostream & os, lpmd::SimulationHistory  *sh) const
{
 // El formato CONFIG de DLPOLY no tiene ningun header especial
 // El formato HISTORY de DLPOLY tiene dos lineas de HEADER.
 if (ftype=="HISTORY")
 {
  if (! Defined("dt")) throw PluginError("dlpoly", "HISTORY format needs the \"dt\" parameter");
  char * buf = new char[512];
  int frame = (*sh).Size();
  int atoms = (*sh)[frame-1].Atoms().Size();
  int records = frame*atoms*(2+level) + frame*4 + 2;
  os << "Generated by LPMD dlpoly plugin\n";
  sprintf(buf, "%10d%10d%10d%10d%10d\n", level, pbkey, frame, atoms, records);
  os << buf;
  *(stepcnt) = 0;
 }
 else os << "Generated by LPMD dlpoly plugin\n";
}

void DlPolyFormat::WriteCell(std::ostream & out, Configuration & con) const
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
#warning SortBySpecies, falta.
 con.SetTag(con, Tag("level"), level);
 char * buf = new char[512];
 if (ftype=="CONFIG")
 {
  sprintf(buf, "%10d%10d\n", level, pbkey);
  out << buf;
  for (int i=0;i<3;++i)
  {
   sprintf(buf, "%20.10f%20.10f%20.10f\n", cell[i][0], cell[i][1], cell[i][2]);
   out << buf; 
  }
  for (long i=0;i<atoms.Size();++i)
  {
   const Atom & at = atoms[i];
   sprintf(buf, "%-5s%13lu\n", at.Symbol().c_str(), i+1);
   out << buf;
   Vector v;
   if (level >= 0)
   {
    v = at.Position();
    for (int q=0;q<3;++q) v[q] = v[q]-0.5*cell[q].Module(); 
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
    out << buf;
   }
   if (level >= 1)
   {
    v = at.Velocity();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
    out << buf;
   }
   if (level >= 2)
   {
    v = at.Acceleration();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
    out << buf;
   }
  } 
  delete [] buf;
 }
 else if(ftype=="HISTORY")
 {
  //Est primera linea en history posee mas informacion.
  //timestep         0       256         0         1    0.001000    0.000000
  //timestep es un string standard; 0 es el current timstep (poner un counter?)
  //256 es el numero de atomos; 0 es el level; 1 es la peridicidad; 0.001000 es el paso del integrador
  //0.000000 es el paso de tiempo real (integrador por current timestep).
  (*stepcnt)++;
  sprintf(buf, "timestep%10ld%10d%10d%10d%12.6g    %12.6g\n", *(stepcnt), int(con.Atoms().Size()), level, pbkey, dt, *(stepcnt)*dt);
  out << buf;
  for (int i=0;i<3;++i)
  {
   sprintf(buf, "%20.10f%20.10f%20.10f\n", cell[i][0], cell[i][1], cell[i][2]);
   out << buf; 
  }
  for (long i=0;i<atoms.Size();++i)
  {
   const Atom & at = atoms[i];
   sprintf(buf, "%-5s%13lu\n", at.Symbol().c_str(), i+1);
   out << buf;
   Vector v;
   if (level >= 0)
   {
    v = at.Position();
    for (int q=0;q<3;++q) v[q] = v[q]-0.5*cell[q].Module(); 
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
    out << buf;
   }
   if (level >= 1)
   {
    v = at.Velocity();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
    out << buf;
   }
   if (level >= 2)
   {
    v = at.Acceleration();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
    out << buf;
   }
  } 
  delete [] buf;
 }
 else 
 {
  throw PluginError("dlpoly","The fkey="+ftype+" value, is not a valid value.");
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new DlPolyFormat(args); }
void destroy(Plugin * m) { delete m; }

