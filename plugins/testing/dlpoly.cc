//
//
//

#include "dlpoly.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>
#include <sstream>
#include <iomanip>

using namespace lpmd;

DlPolyFormat::DlPolyFormat(std::string args): Module("dlpoly"), dt(0.0)
{
 AssignParameter("each", "1");
 AssignParameter("level", "0");
 AssignParameter("periodicity", "1");
 AssignParameter("replacecell", "false");
 AssignParameter("ftype","CONFIG");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 pbkey = GetInteger("periodicity");
 rcell = GetBool("replacecell");
 ftype = GetString("ftype");
 if (Defined("dt")) dt = GetDouble("dt");
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
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso dlpoly.                                     \n";
 std::cout << "      file          : Especifica el archivo que posee el formato dlpoly.       \n";
 std::cout << "      level         : Se especifica el nivel del formato del archivo, estos son\n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << "      periodicity   : Flag periodic boundary key del formato CONFIG            \n";
 std::cout << "      ftype         : File type, puede ser CONFIG o HISTORY.                   \n";
// std::cout << "      megatm        : Number of atoms in simulation cell in last frame. Only   \n";
// std::cout << "                      request for ftype=HISTORY in output mode.                \n";
// std::cout << "      frame         : Number configuration frames in file. Only request for    \n";
// std::cout << "                      ftype=HISTORY in output mode.                            \n";
// std::cout << "      records       : Number of records in file. Only request for ftype=HISTORY\n";
// std::cout << "                      in output mode.                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=dlpoly file=CONFIG level=0                                       \n";
 std::cout << " output module=dlpoly file=CONFIG.output level=1 periodicity=2 each=5          \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato DLPOLY, en el  \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string DlPolyFormat::Keywords() const
{
 return "file each level periodicity replacecell ftype";
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
bool DlPolyFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 sc.MetaData().AssignParameter("level",ToString<int>(level));
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
   if (GetString("replacecell") == "true") sc.SetVector(i, cv[i]);
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
    if (i == 0) this_atom.SetPos(Vector(x+0.5*cv[0].Mod(), y+0.5*cv[1].Mod(), z+0.5*cv[2].Mod()));
    if (i == 1) this_atom.SetVel(Vector(x, y, z));
    if (i == 2) this_atom.SetAccel(Vector(x, y, z));
   }
   sc.AppendAtom(this_atom);
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
   if(GetString("replacecell") == "true") sc.SetVector(j,cv[j]);
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
    if (i == 0) this_atom.SetPos(Vector(x+0.5*cv[0].Mod(), y+0.5*cv[1].Mod(), z+0.5*cv[2].Mod()));
    if (i == 1) this_atom.SetVel(Vector(x, y, z));
    if (i == 2) this_atom.SetAccel(Vector(x, y, z));
   }
   sc.AppendAtom(this_atom);
  }
 }
 else
 {
  throw PluginError("dlpoly", "The fkey="+ftype+" is not a valid fkey value.");
 }
 return true;
}

void DlPolyFormat::WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell>  *cell) const
{
 // El formato CONFIG de DLPOLY no tiene ningun header especial
 // El formato HISTORY de DLPOLY tiene dos lineas de HEADER.
 if (ftype=="HISTORY")
 {
  if (! Defined("dt")) throw PluginError("dlpoly", "HISTORY format needs the \"dt\" parameter");
  char * buf = new char[512];
  int frame = (*cell).size();
  int atoms = (*cell)[frame-1].Size();
  int records = frame*atoms*(2+level) + frame*4 + 2;
  os << "Generated by LPMD dlpoly plugin\n";
  sprintf(buf, "%10d%10d%10d%10d%10d\n", level, pbkey, frame, atoms, records);
  os << buf;
  *(stepcnt) = 0;
 }
 else os << "Generated by LPMD dlpoly plugin\n";
}

void DlPolyFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 sc.SortBySpecies();
 sc.MetaData().AssignParameter("level",ToString<int>(level));
 char * buf = new char[512];
 if (ftype=="CONFIG")
 {
  sprintf(buf, "%10d%10d\n", level, pbkey);
  out << buf;
  for (int i=0;i<3;++i)
  {
   sprintf(buf, "%20.10f%20.10f%20.10f\n", sc.GetVector(i).Get(0), sc.GetVector(i).Get(1), sc.GetVector(i).Get(2));
   out << buf; 
  }
  for (long i=0;i<sc.Size();++i)
  {
   const Atom & at = sc[i];
   sprintf(buf, "%-5s%13lu\n", at.Symb().c_str(), i+1);
   out << buf;
   Vector v;
   if (level >= 0)
   {
    v = at.Position();
    for (int q=0;q<3;++q) v.Set(q, v.Get(q)-0.5*sc.GetVector(q).Mod()); 
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v.Get(0), v.Get(1), v.Get(2));
    out << buf;
   }
   if (level >= 1)
   {
    v = at.Velocity();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v.Get(0), v.Get(1), v.Get(2));
    out << buf;
   }
   if (level >= 2)
   {
    v = at.Acceleration();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v.Get(0), v.Get(1), v.Get(2));
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
  sprintf(buf, "timestep%10ld%10d%10d%10d%12.6g    %12.6g\n", *(stepcnt), sc.Size(), level, pbkey, dt, *(stepcnt)*dt);
  out << buf;
  for (int i=0;i<3;++i)
  {
   sprintf(buf, "%20.10f%20.10f%20.10f\n", sc.GetVector(i).Get(0), sc.GetVector(i).Get(1), sc.GetVector(i).Get(2));
   out << buf; 
  }
  for (long i=0;i<sc.Size();++i)
  {
   const Atom & at = sc[i];
   sprintf(buf, "%-5s%13lu\n", at.Symb().c_str(), i+1);
   out << buf;
   Vector v;
   if (level >= 0)
   {
    v = at.Position();
    for (int q=0;q<3;++q) v.Set(q, v.Get(q)-0.5*sc.GetVector(q).Mod()); 
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v.Get(0), v.Get(1), v.Get(2));
    out << buf;
   }
   if (level >= 1)
   {
    v = at.Velocity();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v.Get(0), v.Get(1), v.Get(2));
    out << buf;
   }
   if (level >= 2)
   {
    v = at.Acceleration();
    sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v.Get(0), v.Get(1), v.Get(2));
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
Module * create(std::string args) { return new DlPolyFormat(args); }
void destroy(Module * m) { delete m; }

