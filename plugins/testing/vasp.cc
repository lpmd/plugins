//
//
//

#include "vasp.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>
#include <sstream>
#include <iomanip>

using namespace lpmd;

VaspFormat::VaspFormat(std::string args): Module("vasp")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("file");
 DefineKeyword("species", "NULL");
 DefineKeyword("each", "1");
 DefineKeyword("level", "0");
 DefineKeyword("type", "Direct");
 AssignParameter("replacecell", "false");
 ProcessArguments(args);
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 speclist = GetString("species");
 satoms = SplitTextLine(speclist,',');
 tp = GetString("type");
 rcell = GetBool("replacecell");
}

VaspFormat::~VaspFormat() { }

void VaspFormat::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la lectura/escritura de archivos en formato  \n";
 std::cout << " POSCAR de VASP, este es un formato con posiciones atomicas.                   \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso vasp.                                       \n";
 std::cout << "      file          : Especifica el archivo que posee el formato POSCAR.       \n";
 std::cout << "      level         : Se especifica el nivel del formato del archivo, estos son\n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << "      species       : Lista las especies (en orden) del fichero POSCAR.        \n";
 std::cout << "                      utilizado en input.                                      \n";
 std::cout << "      type          : Tipo de red Directa/Cartesian para POSCAR.               \n";
 std::cout << "                      utilizado en output.                                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=vasp file=POSCAR level=0                                         \n";
 std::cout << " output module=vasp file=POSCAR.output level=1 each=5                          \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato POSCAR, en el  \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
}

void VaspFormat::ReadHeader(std::istream & is) const { }

// 
// Lee una configuracion desde un archivo CONFIG 
//
bool VaspFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 if(speclist=="NULL") throw PluginError("vasp", "Error, you must pass a species list.");
 double scale=1.0e0;
 std::string tmp;
 Vector cv[3];
 double x, y, z;
 getline(is, tmp);           // lee la linea de titulo
 if (is.eof()) return false; // no hay mas configuraciones que leer
 //lee factor de escala
 getline(is, tmp);
 std::istringstream ost(tmp);
 ost >> scale;
 //lee vectores base cv[]
 for (int i=0;i<3;++i)
 {
  getline(is, tmp);
  std::istringstream vst(tmp);
  vst >> x >> y >> z;
  cv[i] = scale*Vector(x, y, z);
  if (GetString("replacecell") == "true") sc.SetVector(i, cv[i]);
 } 
 //lee y chequea especies atomicas.
 getline(is,tmp);
 RemoveUnnecessarySpaces(tmp);
 std::vector<std::string> numesp = SplitTextLine(tmp,' ');

 //lee si la red es directa o cartesiana, el caso de selective dynamics lo ignora.
 getline(is,tmp);
 RemoveUnnecessarySpaces(tmp);
 std::string tipo="Direct";
 if(tmp[0]=='s' || tmp[0]=='S') getline(is,tmp);

 if(tmp[0]=='c' || tmp[0]=='C') 
 {
  tipo="Cartesian";
 }
 else if(tmp[0]=='d' || tmp[0]=='D')
 {
  tipo="Direct";
 }
 else ShowWarning("plugin vasp", "The type of cell couldn't be read correctly, type=\"Direct\" was assumed.");

 if (numesp.size()!=satoms.size()) { throw PluginError("vasp", "Error, the species number and file POSCAR are different."); }
 else 
 {
  long int S=numesp.size();
  for(int i=0;i<S;++i)
  {
   std::string symbol=satoms[i];
   int ns=(int)atof(numesp[i].c_str());
   for(int j=0;j<ns;++j)
   {
    getline(is,tmp);
    std::istringstream vst(tmp);
    vst >> x >> y >> z;
    Vector vtmp = Vector(x, y, z);
    Atom * this_atom = new Atom(ElemNum(symbol));
    if (tipo=="Cartesian")
    {
     this_atom->SetPos(scale*vtmp);
    }
    else if (tipo=="Direct")
    {
     vtmp = (vtmp[0]*cv[0],vtmp[1]*cv[1],vtmp[2]*cv[2]);
     this_atom->SetPos(vtmp);
    }
    else
    {
     throw PluginError("vasp","Unexpected error setting atomic position, check 'type'.");
    }
    sc.Create(this_atom);
   }
  }
 }
 return true;
}

void VaspFormat::WriteHeader(std::ostream & os, std::vector<SimulationCell> * cells) const
{
}

void VaspFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 out << "Generated by LPMD vasp plugin\n";
 out << "  1  " << '\n'; 
 for (int i=0;i<3;++i)
 {
  out.setf(std::ios::left);
  out.setf(std::ios::fixed);
  out << " " << std::setw(8) << std::setprecision(8) << sc.GetVector(i)[0]; 
  out << " " << std::setw(8) << std::setprecision(8) << sc.GetVector(i)[1]; 
  out << " " << std::setw(8) << std::setprecision(8) << sc.GetVector(i)[2]; 
  out << '\n';
 }
 std::list<std::string> symesp=sc.SpeciesList();
 std::vector<int> nspec = sc.NSpecies();
 for (unsigned long i=0;i<nspec.size();i++) 
 {
  out << " " << nspec[i] << " ";
 }
 out << '\n';
 out << tp << '\n';
 for(std::list<std::string>::const_iterator it=symesp.begin();it!=symesp.end();++it)
 {
  for (unsigned long int i=0;i<sc.size();++i)
  {
   const Atom & at = sc[i];
   if(at.Symb()==*it)
   {
    if(tp=="Cartesian")
    {
     out << at.Position() << '\n';
    }
    else if(tp=="Direct")
    {
     Vector tmp = at.Position();
     double lx = (sc.GetVector(0)).Module();
     double ly = (sc.GetVector(1)).Module();
     double lz = (sc.GetVector(2)).Module();
     tmp = Vector(tmp[0]/lx,tmp[1]/ly,tmp[2]/lz);
     out << tmp << '\n';
    }
   }
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VaspFormat(args); }
void destroy(Module * m) { delete m; }

