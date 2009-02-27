//
//
//

#include "pdb.h"

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>
#include <algorithm>
#include <cctype>

using namespace lpmd;

PDBFormat::PDBFormat(std::string args): Module("pdb")
{
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("file");
 DefineKeyword("each", "1");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 writefile = GetString("file");
 interval = GetInteger("each");
}

PDBFormat::~PDBFormat() { }

void PDBFormat::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la escritura de archivos en formato PDB     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso pdb.                                        \n";
 std::cout << "      file          : Especifica el archivo que posee el formato pdb.          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " output module=pdb file=outputfile.pdb each=5                                \n\n";
 std::cout << "      De esta forma podemos escribir archivos en formato pdb.                  \n";
}

void PDBFormat::WriteHeader(std::ostream & os, std::vector<SimulationCell> * cells) const
{
 // PDB no tiene ningun header especial
}

void PDBFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 out << "COMPND   From lpmd\n";
 out << "AUTHOR   pdb plugin for lpmd\n";
 for (unsigned long int i=0;i<sc.size();++i)
 {
  const Vector & pos = sc[i].Position();
  out << std::setw(6) << "HETATM"; //1-6
  out << std::setw(5) << (i+1); //7-11
  std::string sy = sc[i].Symb();
  std::transform(sy.begin(), sy.end(), sy.begin(), (int(*) (int)) std::toupper);
  out << std::setw(1) << " ";
  out << std::setw(4) <<std::left<< sy ; //13-16
  out << std::setw(1) << " "; //17
  out << std::setw(3) <<std::left<< "LIG"; //18-20
  out << std::right;
  out << std::setw(1) << " ";
  out << std::setw(1) << " "; //22
  out << std::setw(4) << (i+1); //23-26
  out << std::setw(1) << " "; //27
  out << std::setw(3) << "  ";
  out << std::setw(8) << std::fixed <<std::setprecision(3)<< pos.GetX(); //31-38
  out << std::setw(8) << std::fixed <<std::setprecision(3)<< pos.GetY(); //39-46
  out << std::setw(8) << std::fixed <<std::setprecision(3)<< pos.GetZ(); //47-54
  out << std::right;
  out << std::setw(7) << "1.00"; //55-60
  out << std::setw(6) << "0.00"; //61-66
  out << std::setw(8) << "        ";
  out << std::setw(4) << " "; //73-76
  out << std::setw(2) << sy ; //77-78
  out << std::setw(2) << " "; //79-80
  out << '\n';
 }
 out << "END\n\n";
}

void PDBFormat::ReadHeader(std::istream & is) const
{
 // PDB no tiene header especial para leer
}

bool PDBFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 long natoms=0;
 std::vector<lpmd::Atom> atomlist;
 if(is.eof()) return false;
 while(!is.eof())
 {
  std::string tmp;
  getline(is,tmp);
  if(tmp[0]!='#')
  {
   if (tmp.compare(0,6,"HETATM")==0)
   {
    natoms++;
    std::vector<std::string> linea = SplitTextLine(tmp,' ');
    std::string symb = linea[2];
    if(symb.length()>1)
    {
     for(unsigned int j=1;j<symb.length();j++) symb[j]=tolower(symb[j]);
    }
    lpmd::Vector pos(atof(linea[5].c_str()),atof(linea[6].c_str()),atof(linea[7].c_str()));
    int elem=ElemNum(symb);
    sc.Create(new Atom(elem,pos));
   }
   if (tmp.compare(0,3,"END")==0)
   {
    return true;
   }
  }
 }
 return false;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new PDBFormat(args); }
void destroy(Module * m) { delete m; }

