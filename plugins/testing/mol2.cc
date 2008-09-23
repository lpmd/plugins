//
//
//

#include "mol2.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>
#include <algorithm>
#include <cctype>

using namespace lpmd;

Mol2Format::Mol2Format(std::string args): Module("mol2")
{
 AssignParameter("each", "1");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 writefile = GetString("file");
 interval = GetInteger("each");
 rcell = GetBool("replacecell");
}

Mol2Format::~Mol2Format() { }

void Mol2Format::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = mol2                                                     \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la escritura de archivos en formato Mol2     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso mol2.                                       \n";
 std::cout << "      file          : Especifica el archivo que posee el formato lpmd.         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " output module=mol2 file=outputfile.mol2 each=5                              \n\n";
 std::cout << "      De esta forma podemos escribir archivos en formato mol2.                 \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Mol2Format::Keywords() const 
{
 return "file each replacecell";
}

void Mol2Format::WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *cell) const
{
 // Mol2 no tiene ningun header especial
}

void Mol2Format::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 out << "@<TRIPOS>MOLECULE\n";
 out << "autogenerated\n";
 out << " " << sc.Size() << " 0 0 0 0\n";
 out << "SMALL\n";
 out << "GASTEIGER\n";
 out << "Energy = 0\n\n";
 out << "@<TRIPOS>ATOM\n";
 for (long i=0;i<sc.Size();++i)
 {
  const Vector & pos = sc[i].Position();
  out << std::setw(8) << (i+1);
  std::string sy = sc[i].Symb();
  std::transform(sy.begin(), sy.end(), sy.begin(), (int(*) (int)) std::toupper);
  out << "  " << sy << "    ";
  out << pos << " " << sc[i].Symb() << "   0  LIG0      0.0000\n";
 }
 out << "@<TRIPOS>BOND\n\n";
}

void Mol2Format::ReadHeader(std::istream & is) const
{
 // Mol2 no tiene header especial para leer
}

bool Mol2Format::ReadCell(std::istream & is, SimulationCell & sc) const
{
 if (GetString("replacecell") == "true") throw PluginError("mol2", "This format does not contain any cell vectors.");
 long natoms=0;
 if(is.eof()) return false;
 while(!is.eof())
 {
  std::string tmp;
  getline(is,tmp);
  std::cerr << tmp << '\n';
  if(tmp[0]!='#')
  {
   if (tmp=="@<TRIPOS>MOLECULE")
   {
    getline(is,tmp);
    getline(is,tmp);
    std::vector<std::string> linea = SplitTextLine(tmp,' ');
    natoms = (long)atof(linea[0].c_str()); 
    std::cerr << " numero de atomos >> " << natoms << '\n';
   }
   if (tmp=="@<TRIPOS>ATOM")
   {
    if (natoms <=0) throw PluginError("mol2","Atoms number in mol2 file was not be read!");
    else
    {
     for (int i=0;i<natoms;i++)
     {
      getline(is,tmp);
      std::vector<std::string> linea = SplitTextLine(tmp,' ');
      std::string symb = linea[1];
      lpmd::Vector pos(atof(linea[2].c_str()),atof(linea[3].c_str()),atof(linea[4].c_str()));
      if(symb.length()>1)
      {
       for(unsigned int j=1;j<symb.length();j++)
       {
	symb[j]=tolower(symb[j]);
       }
      }
      int elem=ElemNum(symb);
      std::cerr << "Anadiendo atomo " << elem << " - " << pos << '\n';
      sc.AppendAtom(Atom(elem,pos));
     }
    }
    return true;
   }
  }
 }
 return false;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Mol2Format(args); }
void destroy(Module * m) { delete m; }

