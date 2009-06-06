//
//
//

#include "dlpoly.h"

#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <sstream>
#include <iomanip>

using namespace lpmd;

DlPolyFormat::DlPolyFormat(std::string args): Plugin("dlpoly", "2.0")
{
 ParamList & param = (*this);
 //
 DefineKeyword("file");
 DefineKeyword("each", "1");
 DefineKeyword("level", "0");
 DefineKeyword("periodicity", "1");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = param["file"];
 interval = int(param["each"]);
 level = int(param["level"]);
 pbkey = int(param["periodicity"]);
 rcell = bool(param["replacecell"]);
}

DlPolyFormat::~DlPolyFormat() { }

void DlPolyFormat::ShowHelp() const
{
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
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=dlpoly file=CONFIG level=0                                       \n";
 std::cout << " output module=dlpoly file=CONFIG.output level=1 periodicity=2 each=5          \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato DLPOLY, en el  \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
}

void DlPolyFormat::ReadHeader(std::istream & is) const
{
 // El formato CONFIG de DLPOLY no tiene ningun header especial
}

// 
// Lee una configuracion desde un archivo CONFIG 
//
bool DlPolyFormat::ReadCell(std::istream & is, Configuration & con) const
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 int fkey, pbk;
 std::string tmp;
 Vector cv[3];
 double x, y, z;
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
  Atom this_atom(symbol);
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
 return true;
}

void DlPolyFormat::WriteHeader(std::ostream & os, SimulationHistory * sh) const
{
 // El formato CONFIG de DLPOLY no tiene ningun header especial
}

void DlPolyFormat::WriteCell(std::ostream & out, Configuration & con) const
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 char * buf = new char[512];
 out << "Generated by LPMD dlpoly plugin\n";
 sprintf(buf, "%10d%10d\n", level, pbkey);
 out << buf;
 for (int i=0;i<3;++i)
 {
  sprintf(buf, "%20.10f%20.10f%20.10f\n", cell[i][0], cell[i][1], cell[i][2]);
  out << buf; 
 }
 for (long int i=0;i<atoms.Size();++i)
 {
  sprintf(buf, "%-5s%13lu\n", atoms[i].Symbol().c_str(), i+1);
  out << buf;
  Vector v;
  if (level >= 0)
  {
   v = atoms[i].Position();
   for (int q=0;q<3;++q) v[q] = v[q]-0.5*cell[q].Module(); 
   sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
   out << buf;
  }
  if (level >= 1)
  {
   v = atoms[i].Velocity();
   sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
   out << buf;
  }
  if (level >= 2)
  {
   v = atoms[i].Acceleration();
   sprintf(buf, "    %12.10g        %12.10g        %12.10g\n", v[0], v[1], v[2]);
   out << buf;
  }
 } 
 delete [] buf;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new DlPolyFormat(args); }
void destroy(Plugin * m) { delete m; }

