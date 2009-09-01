//
//
//

#include "rawbinary.h"

#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <lpmd/atom.h>

using namespace lpmd;

RawBinFormat::RawBinFormat(std::string args): Plugin("rawbinary", "2.0")
{
 ParamList & param = (*this);
 AssignParameter("level", "0");
 AssignParameter("each", "1");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = param["file"];
 interval = int(param["each"]);
 level = int(param["level"]);
 rcell = bool(param["replacecell"]);
}

RawBinFormat::~RawBinFormat() { }

void RawBinFormat::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la lectura/escritura de archivos en formato  \n";
 std::cout << " binario.                                                                      \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso rawbinary.                                  \n";
 std::cout << "      file          : Especifica el archivo que posee el formato lpmd.         \n";
 std::cout << "      level         : Se especifica el nivel del formato de rawbinary, estos son \n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=rawbinary file=inputfile.raw level=0                             \n";
 std::cout << " output module=rawbinary file=outputfile.raw level=1 each=5                    \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato rawbinary, en  \n";
 std::cout << "      el caso de la salida, es necesaria la opcion each.                       \n";
}

void RawBinFormat::ReadHeader(std::istream & is) const
{
 assert(&is != 0); //icc 869
 // El formato RawBinary no tiene ningun header especial
}

// 
// Lee una configuracion desde un archivo RawBinary
//
bool RawBinFormat::ReadCell(std::istream & is, Configuration & con) const
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 int level = int(Parameter(con.GetTag(con, Tag("level"))));
 if ((*this)["replacecell"] == "true") throw PluginError("rawbinary", "This format does not contain any cell vectors.");
 long int s = 0;
 int lvl = level;
 is.read((char *)(&s), sizeof(long int));
 if (is.eof()) return false;
 is.read((char *)(&lvl), sizeof(int));
 for (long int i=0;i<s;++i)
 {
  double p;
  int symbol = 0; 
  is.read((char *)(&symbol), sizeof(int));
  if (is.eof()) throw PluginError("rawbinary", "Unexpected end of file on reading");
  Vector pos, vel, acc;
  for (int q=0;q<3;++q)
  {
   is.read((char *)(&p), sizeof(double));
   pos[q] =  p;
  }
  if (lvl > 0)
  {
   for (int q=0;q<3;++q)
   {
    is.read((char *)(&p), sizeof(double));
    vel[q] = p;
   }
  }
  if (lvl > 1)
  {
   for (int q=0;q<3;++q)
   {
    is.read((char *)(&p), sizeof(double));
    acc[q] = p;
   }
  }
  atoms.Append(Atom(ElemSym[symbol], pos, vel, acc));
 }
 return true;
}

void RawBinFormat::WriteHeader(std::ostream & os, lpmd::SimulationHistory * sh) const
{
 assert (&os != 0); //icc 869
 assert (sh >=(void *)NULL); //icc 869
 // El formato RawBinary no tiene ningun header especial
}

void RawBinFormat::WriteCell(std::ostream & out, lpmd::Configuration & con) const
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 int level = int(Parameter(con.GetTag(con, Tag("level"))));
 long int totsize = 0, expsize;
 long int s = atoms.Size();
 expsize = sizeof(long int)+sizeof(int)+s*(sizeof(int)+3*sizeof(double)*(level+1));
 char * buffer = new char[expsize];
 memcpy((void *)(&buffer[totsize]), (void *)(&s), sizeof(long int));
 totsize += sizeof(long int);
 int lvl = level;
 memcpy((void *)(&buffer[totsize]), (void *)(&lvl), sizeof(int));
 totsize += sizeof(int);
 for (long int i=0;i<s;++i)
 {
  double p;
  int symbol = atoms[i].Z();
  memcpy((void *)(&buffer[totsize]), (void *)(&symbol), sizeof(int));
  totsize += sizeof(int);
  Vector pos = atoms[i].Position();
  for (int q=0;q<3;++q)
  {
   p = pos[q];
   memcpy((void *)(&buffer[totsize]), (void *)(&p), sizeof(double));
   totsize += sizeof(double);
  }
  if (level > 0) 
  {
   Vector vel = atoms[i].Velocity();
   for (int q=0;q<3;++q)
   {
    p = vel[q]; 
    memcpy((void *)(&buffer[totsize]), (void *)(&p), sizeof(double));
    totsize += sizeof(double);
   }   
  }
  if (level > 1)
  {
   Vector acc = atoms[i].Acceleration();
   for (int q=0;q<3;++q)
   {
    p = acc[q]; 
    memcpy((void *)(&buffer[totsize]), (void *)(&p), sizeof(double));
    totsize += sizeof(double);
   }   
  } 
 }
 if (expsize != totsize)
 {
  delete [] buffer;
  throw PluginError("rawbinary", "Unexpected byte mismatch on writing");
 }
 else
 {
  out.write((const char *)(buffer), expsize);
  delete [] buffer;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new RawBinFormat(args); }
void destroy(Plugin * m) { delete m; }

