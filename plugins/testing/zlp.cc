//
//
//

#include "zlp.h"

#include <arpa/inet.h>
#include <zlib.h>

#include <lpmd/util.h>
#include <lpmd/simulationcell.h>

#include <sstream>

#define ZLP_VERSION_MAJOR 1
#define ZLP_VERSION_MINOR 0
#define ZLP_VERSION_REV   0

#define ZLP_NONE 0
#define ZLP_READ 1
#define ZLP_WRITE 2

using namespace lpmd;

ZLPFormat::ZLPFormat(std::string args): Module("zlp")
{
 AssignParameter("each", "1");
 AssignParameter("level", "0");
 AssignParameter("blocksize", "1024");
 AssignParameter("compression", "6");
 AssignParameter("replacecell", "false");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 readfile = writefile = GetString("file");
 interval = GetInteger("each");
 level = GetInteger("level");
 blocksize = GetInteger("blocksize");
 complev = GetInteger("compression");
 rcell = GetBool("replacecell");
 // inicializa la estructura z_stream
 zstr = (void *)(new z_stream);
 lastop = new int(ZLP_NONE);
 inbuf = new unsigned char[blocksize];
 outbuf = new unsigned char[blocksize];
}

ZLPFormat::~ZLPFormat() 
{
 if ((*lastop) == ZLP_WRITE) deflate((z_stream *)(zstr), Z_FINISH);
 if ((*lastop) == ZLP_READ) inflateEnd((z_stream *)(zstr));
 if ((*lastop) == ZLP_WRITE) deflateEnd((z_stream *)(zstr));
 delete lastop;
 delete [] inbuf;
 delete [] outbuf;
 delete (z_stream *)(zstr); 
}

void ZLPFormat::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = zlp                                                      \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para la lectura/escritura de archivos en formato  \n";
 std::cout << " zlp, este es un formato comprimido con posiciones escaladas y propio de lpmd. \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      module        : En la opcion input, es necesario especificar el formato  \n";
 std::cout << "                      en este caso zlp.                                        \n";
 std::cout << "      file          : Especifica el archivo que posee el formato zlp.          \n";
 std::cout << "      level         : Se especifica el nivel del formato de zlp, estos son     \n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-ace.                       \n";
 std::cout << "      blocksize     : Especifica el tamanyo del buffer interno de compresion.  \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " input module=zlp file=inputfile.zlp level=0                                   \n";
 std::cout << " output module=zlp file=outputfile.zlp level=1 each=5                          \n\n";
 std::cout << "      De esta forma podemos leer o escribir archivos en formato zlp, en el     \n";
 std::cout << " en el caso de la salida, es necesaria la opcion each.                         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string ZLPFormat::Keywords() const 
{
 return "file each level blocksize replacecell";
}

void ZLPFormat::ReadHeader(std::istream & is) const
{ 
 char hdr[10];
 is.read(hdr, 10);
 if (((hdr[0] != 'Z') || (hdr[2] != 'L')) || (hdr[4] != 'P')) throw PluginError("zlp", "Wrong header");
 if ((hdr[1] != '4') || (hdr[3] != '2')) throw PluginError("zlp", "Wrong header");
 unsigned short int v0 = hdr[5];   // major version number
 unsigned short int v1 = hdr[6];   // minor version number
 unsigned short int v2 = hdr[7];   // revision 
 unsigned short int bf0 = hdr[8];  // reservado para 8-bit flag
 unsigned short int bf1 = hdr[9];  // reservado para 8-bit flag
 std::cerr << "[zlp] version number from file is " << v0 << "." << v1 << "." << v2 << '\n';
 std::cerr << "[zlp] format flags are " << bf0 << " and " << bf1 << '\n';
 z_stream & stream = *((z_stream *)(zstr));
 stream.zalloc = Z_NULL;
 stream.zfree = Z_NULL;
 stream.opaque = Z_NULL;
 if (inflateInit(&stream) != Z_OK) throw PluginError("zlp", "Decompression failed");
}

// 
// Reads a configuration from a ZLP file 
//
bool ZLPFormat::ReadCell(std::istream & is, SimulationCell & sc) const
{
 sc.MetaData().AssignParameter("level",ToString<int>(level));
 *lastop = ZLP_READ;
 z_stream & stream = *((z_stream *)(zstr));
 std::istringstream bufstr(std::istringstream::in);
 long int ts = 0;                                // total de bytes leidos hasta ahora
 std::string * istr = new std::string;
 unsigned char foo;
 unsigned long int cbufs;
 is.read((char *)&foo, 1);
 if (is.eof()) return false;
 is.read((char *)&cbufs, int(foo));
 cbufs = ntohl(cbufs);
 // std::cerr << "DEBUG This configuration has " << cbufs << " compressed bytes" << '\n';
 while (1)
 {
  long int rem = cbufs - ts;
  if (rem == 0) break;
  long int chunk = blocksize;
  if (chunk > rem) chunk = rem;
  is.read((char *)inbuf, chunk);
  ts += is.gcount();
  //std::cerr << "DEBUG read " << is.gcount() << " bytes from ZLP file" << '\n';
  stream.avail_in = is.gcount();
  stream.next_in = (unsigned char *)(inbuf);
  //std::cerr << "DEBUG decompressing..." << '\n';
  do
  {
   stream.avail_out = blocksize; 
   stream.next_out = (unsigned char *)(outbuf);
   int f = (is.eof() ? Z_FINISH : Z_SYNC_FLUSH);
   inflate(&stream, f);                           // no chequea el estado aun 
   int have = blocksize-stream.avail_out;
   //std::cerr << "DEBUG read chunk of " << have << " uncompressed bytes" << '\n';
   for (int q=0;q<have;++q) (*istr) += (char)(outbuf[q]);
  } while (stream.avail_out == 0);
 }
 std::cerr << "-> ZLP decompression: expanded " << cbufs << " bytes to " << istr->size() << " bytes.\n";
 //
 //
 // 
 std::istringstream ibufstr(*istr);
 int lvl;
 long int natoms;
 ibufstr >> lvl;
 ibufstr >> natoms;
 //std::cerr << "DEBUG Number of atoms = " << natoms << '\n';
 for (int j=0;j<3;++j)
 { 
  Vector v;
  double vq[3];
  for (int i=0;i<3;++i) 
  {
   ibufstr >> vq[i];
   v.Set(i, vq[i]);
  }
  if (GetString("replacecell") == "true") sc.SetVector(j, v);
 }
 for (long int i=0;i<natoms;++i) 
 {
  Atom at;
  std::string sym;
  ibufstr >> sym;
  at.SetSpc(ElemNum(sym));
  Vector fpos, vel, acc;
  double vq[3];
  for (int q=0;q<3;++q) 
  {
   ibufstr >> vq[q];
   fpos.Set(q, vq[q]);
  }
  sc.AppendAtom(at);
  sc.SetFracPosition(i, fpos);
  if (lvl > 0)
  {
   for (int q=0;q<3;++q) 
   {
    ibufstr >> vq[q];
    vel.Set(q, vq[q]);
   }
   sc.SetVelocity(i, vel);
  }
  if (lvl > 1) 
  {
   for (int q=0;q<3;++q) 
   {
    ibufstr >> vq[q];
    acc.Set(q, vq[q]);
   }
   sc.SetAcceleration(i, acc);
  }
 }
 delete istr;
 return true;
}

void ZLPFormat::WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *cell) const
{
 sc.MetaData().AssignParameter("level",ToString<int>(level));
 char hdr[10] = {'Z', '4', 'L', '2', 'P', 1, 0, 0, 0, 0};
 // Z4L2P es la firma que marca el archivo como formato ZLP
 // los tres bytes siguientes son los numeros de version del formato
 // los dos ultimos bytes son reservados para 16 flags booleanos
 // 
 // aqui se deberian setear los 16 flags
 // 
 // se escriben los 10 bytes del header
 os.write(hdr, 10);
 z_stream & stream = *((z_stream *)(zstr));
 stream.zalloc = Z_NULL;
 stream.zfree = Z_NULL;
 stream.opaque = Z_NULL;
 if (deflateInit(&stream, complev) != Z_OK) throw PluginError("zlp", "Compression failed");
}

//
// Writes a configuration to a ZLP file
//
void ZLPFormat::WriteCell(std::ostream & out, SimulationCell & sc) const
{
 *lastop = ZLP_WRITE;
 z_stream & stream = *((z_stream *)(zstr));
 std::ostringstream * ostr = new std::ostringstream();
 std::ostringstream & obufstr = *ostr;
 // Escribe la configuracion completa en bufstr 
 // Esto gasta memoria enormemente pero por ahora deberia bastar
 obufstr << level << '\n';
 obufstr << sc.Size() << '\n';
 for (int j=0;j<3;++j) 
   for (int i=0;i<3;++i) { obufstr.precision(15); obufstr << std::fixed << sc.GetVector(j).Get(i) << '\n'; }
 for (long int i=0;i<sc.Size();++i) 
 {
  obufstr << sc.GetAtom(i).Symb() << '\n';
  for (int q=0;q<3;++q) { obufstr.precision(15); obufstr << std::fixed << sc.FracPosition(i).Get(q) << '\n'; }
  if (level > 0) 
     for (int q=0;q<3;++q) { obufstr.precision(15); obufstr << std::fixed << sc.GetAtom(i).Velocity().Get(q) << '\n'; }
  if (level > 1) 
     for (int q=0;q<3;++q) { obufstr.precision(15); obufstr << std::fixed << sc.GetAtom(i).Acceleration().Get(q) << '\n'; }
 }
 std::istringstream * istr = new std::istringstream(obufstr.str());
 delete ostr;
 std::istringstream & ibufstr = *istr;
 //
 //
 //
 long int ucsize = ibufstr.str().size();
 std::string cbuf;
 while (1)
 {
  ibufstr.read((char *)inbuf, blocksize);
  //std::cerr << "DEBUG writing " << ibufstr.gcount() << " uncompressed bytes to ZLP file" << '\n';
  if (ibufstr.gcount() == 0) break;
  stream.avail_in = ibufstr.gcount();
  stream.next_in = inbuf;
  do 
  {
   stream.avail_out = blocksize;
   stream.next_out = outbuf;
   if (deflate(&stream, Z_SYNC_FLUSH) == Z_STREAM_ERROR) throw PluginError("zlp", "Compression failed");
   int have = blocksize-stream.avail_out;
   for (int q=0;q<have;++q) cbuf += (char)(outbuf[q]);
   //std::cerr << "DEBUG writing " << have << " compressed bytes to ZLP file" << '\n';
  } while (stream.avail_out == 0);
 }
 //std::cerr << "DEBUG compressed data for this configuration is " << cbuf.size() << " bytes\n";
 std::cerr << "-> ZLP compression: packed " << ucsize << " bytes into " << cbuf.size() << " bytes.\n";
 unsigned char foo = sizeof(unsigned long int);
 unsigned long int cbufs = htonl(cbuf.size());
 out.write((char *)&foo, 1);
 out.write((char *)&cbufs, int(foo));
 out.write(cbuf.c_str(), cbuf.size());
 out.flush();
 // 
 //
 delete istr;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ZLPFormat(args); }
void destroy(Module * m) { delete m; }

