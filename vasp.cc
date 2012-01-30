//
//
//

#include "vasp.h"

#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <sstream>
#include <iomanip>

using namespace lpmd;

VaspFormat::VaspFormat(std::string args): Plugin("vasp", "3.1")
{
 ParamList & params = (*this);
 //
 DefineKeyword("file");
 DefineKeyword("species", "NULL");
 DefineKeyword("numatoms", "NULL");
 DefineKeyword("each", "1");
 DefineKeyword("level", "0");
 DefineKeyword("type", "Direct");
 DefineKeyword("ftype", "poscar");
 AssignParameter("replacecell", "false");
 ProcessArguments(args);
 readfile = writefile = params["file"];
 interval = int(params["each"]);
 level = int(params["level"]);
 speclist = params["species"];
 numatoms = params["numatoms"];
 satoms = StringSplit(speclist,',');
 numesp = StringSplit(numatoms,',');
 tp = params["type"];
 ftype = params["ftype"];
 rcell = params["replacecell"];
 first=true;
}

VaspFormat::~VaspFormat() { }

void VaspFormat::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = vasp                                                     \n";
 std::cout << " Problems Report to = admin@lpmd.cl                                            \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin is used to read/write VASP's atomic configurations files     \n";
 std::cout << "      (POSCAR files).                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      file          : Input/output file that contains the atomic configurations\n";
 std::cout << "                      in POSCAR format.                                        \n";
 std::cout << "      species       : Atomic species list (in order) of the POSCAR file.       \n";
 std::cout << "      numatoms      : Number of atoms per species.                             \n";
 std::cout << "      each          : Determines how often (each how many time-steps) the      \n";
 std::cout << "                      input/output file must be read/written.                  \n";
 std::cout << "      level         : Determines the file format level (0/1/2).                \n";
 std::cout << "                      0/1/2 <-> pos/pos-vel/pos-vel-acel.                      \n";
 std::cout << "      type          : Positions type (Direct / Cartesian) for the              \n";
 std::cout << "                      POSCAR/CONTCAR file.                                     \n";
 std::cout << "      ftype         : File type that is going to be read/written               \n";
 std::cout << "                      (poscar/contcar/xdatcar):                                \n";
 std::cout << "                      poscar: To read/write POSCAR files. They contain the     \n";
 std::cout << "                              positions of the ions and their velocities       \n";
 std::cout << "                              (if there are any).                              \n";
 std::cout << "                      contcar: To read/write CONTCAR files. They contain the   \n";
 std::cout << "                               same as the POSCAR file, but an extra block     \n";
 std::cout << "                               containing the predictor-corrector coordinates  \n";
 std::cout << "                               needed as an input for the next MD-job in VASP. \n";
 std::cout << "                               LPMD will skip this block when reading, and will\n";
 std::cout << "                               not print it when writing.                      \n";
 std::cout << "                      xdatcar: To read/write XDATCAR files. They contain the   \n";
 std::cout << "                               position of the ions for every ionic-step. Be   \n";
 std::cout << "                               aware that XDATCAR files DO NOT contain cell    \n";
 std::cout << "                               vectors, even when they change.                 \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example           >>                                                          \n";
 std::cout << " #Loading the plugin :                                                         \n";
 std::cout << " input module=vasp file=POSCAR species=Si level=0                              \n";
 std::cout << " output module=vasp file=XDATCAR-01.output species=Si,O ftype=xdatcar numatoms=16,32\n\n";
 std::cout << "      This plugin is used to read and write configurations in VASP's formats   \n";
 std::cout << "      (POSCAR, CONTCAR and XDATCAR). The file extension (.output) is irrelevant,\n";
 std::cout << "      what matters is the module loaded (module=vasp).                         \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

void VaspFormat::ReadHeader(std::istream & is) const {assert(&is != 0);}//icc 869

// 
// Reads file from POSCAR/CONTCAR/XDATCAR
//
bool VaspFormat::ReadCell(std::istream & is, Configuration & con) const
{
 BasicParticleSet & part = con.Atoms();
 BasicCell & cell = con.Cell();
 assert(part.Size()==0);
 if(speclist=="NULL") throw PluginError("vasp", "Error, you must pass a species list.");
 double scale=1.0e0;
 std::string tmp;
 Vector cv[3];
 double x, y, z;
 if (first) getline(is,tmp);  // Read title/atoms for POSCAR/XDATCAR.
 if (is.eof()) return false; // no more configurations to read
 
 // Read Cell Vectors for POSCAR/CONTCAR files
 if (ftype=="poscar" || ftype=="contcar")
 {
  getline(is, tmp);           // read scale factor
  std::istringstream ost(tmp);
  ost >> scale;
  //Read cell vectors, cv[]
  for (int i=0;i<3;++i)
  {
   getline(is, tmp);
   std::istringstream vst(tmp);
   vst >> x >> y >> z;
   cv[i] = scale*Vector(x, y, z);
   if ((*this)["replacecell"] == "true") cell[i] = cv[i];
  }
 }

 // For POSCAR/CONTCAR: Read and check out the amount of atoms for each species
 if (ftype=="poscar" || ftype=="contcar")  getline(is, tmp);
 RemoveUnnecessarySpaces(tmp);
 if (first && numatoms=="NULL") numesp = StringSplit(tmp,' ');

 // Reads the type of atomic positions (Direct/Cartesian). In the case of 'selective dynamics', ignore it.
 std::string tipo="Direct";
 if (ftype=="poscar" || ftype=="contcar")
 {
  getline(is,tmp);
  RemoveUnnecessarySpaces(tmp);
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
 }
 else if (ftype=="xdatcar" && first)
 {
  for (int i=0; i<5; ++i) {getline(is,tmp); first=false;}// skip unnecesary lines
 }

 // Check that the amount of atoms match.
 bool aff = numesp.Size()!=satoms.Size(); // For poscar and contcar, check that N("Si,O")=N("16,32")=2, for example.
 if (ftype=="xdatcar") aff=false;

 if (aff)
 {
  throw PluginError("vasp", "Error, the number of species given does not match with the number of species of the POSCAR file.");
 }
 else 
 {
  long int S=numesp.Size();
  if (ftype=="xdatcar" && numatoms=="NULL") S=1; // Assume one species if none is given.
  for(int i=0;i<S;++i)
  {
   std::string symbol=satoms[i];
   int ns=(int)atof(numesp[i].c_str());
   for(int j=0;j<ns;++j)
   {
    getline(is,tmp);
    RemoveUnnecessarySpaces(tmp);
    std::istringstream vst(tmp);
    vst >> x >> y >> z;
    Vector vtmp = Vector(x, y, z);
    Atom this_atom(symbol);
    if (tipo=="Cartesian")
    {
     this_atom.Position() = scale*vtmp;
    }
    else if (tipo=="Direct")
    {
     vtmp = cell.Cartesian(vtmp); // Antes era: vtmp = vtmp[0]*cv[0]+vtmp[1]*cv[1]+vtmp[2]*cv[2];
     this_atom.Position() = vtmp;
    }
    else
    {
     throw PluginError("vasp","Unexpected error setting atomic position, check 'type'.");
    }
    part.Append(this_atom);
   }
  }


  // Read velocities, if there are.
  if (ftype=="poscar" || ftype=="contcar")
  {
   getline(is,tmp);
   getline(is,tmp);
   RemoveUnnecessarySpaces(tmp);

   if (is.eof()) return false; // no more configurations to read
   while (tmp=="" && !is.eof()) {getline(is,tmp); RemoveUnnecessarySpaces(tmp);}
   
   if (tmp!="")
   {
    int n=0;
    for (int i=0; i<numesp.Size(); ++i) n+=(int)atof(numesp[i].c_str());
    for(int j=0;j<n;++j)
    {
     if (is.eof())
     {
      throw PluginError("vasp", "Error, the number of velocities does not match with the number of atoms.");
      return false;
     }
     std::istringstream vst(tmp);
     vst >> x >> y >> z;
     Vector vtmp = Vector(x, y, z);
     part[j].Velocity() = vtmp;
     getline(is,tmp);
    }
    // Skip everything below the velocities
    while (!is.eof()) getline(is,tmp); 
   }
  }
  else if (ftype=="xdatcar")
  {
   getline(is,tmp);  // Read the blank space below the block and throw an error if it's not blank
   if (!first && StringSplit(tmp,' ').Size()>2) throw PluginError("vasp", "Error, numatoms does not match with the actual number of atoms in the file.");
  }
 }
 con.SetTag(con, Tag("level"), 1);
 return true;
}

void VaspFormat::WriteHeader(std::ostream & os, SimulationHistory * sh) const
{
 assert(&os != 0); //icc 869
 assert(sh >= (void *)NULL); //icc 869
}

void VaspFormat::WriteCell(std::ostream & out, Configuration & con) const
{
 BasicParticleSet & part = con.Atoms();
 BasicCell & cell = con.Cell();

 out << "Generated by LPMD's vasp plugin\n";
 out << "   1.0" << '\n'; 
 for (int i=0;i<3;++i)
 {
  out.setf(std::ios::left);
  out.setf(std::ios::fixed);
  out << "     " << std::setw(8) << std::setprecision(8) << cell[i][0]; 
  out << "     " << std::setw(8) << std::setprecision(8) << cell[i][1]; 
  out << "     " << std::setw(8) << std::setprecision(8) << cell[i][2]; 
  out << '\n';
 }
 Array<int> esp=part.Elements();
 Array<Array <int> > list;

 for (long i=0;i<esp.Size();++i)
 {
  list.Append(part.WithZ(esp[i]));
 }

 for (long i=0;i<list.Size();++i) 
 {
  out << " " << list[i].Size() << " ";
 }
 out << '\n';
 out << tp << '\n';
 for(long i=0 ; i < list.Size() ; ++i)
 {
  for (long int j=0;j<list[i].Size();++j)
  {
   if(tp=="Cartesian")
   {
    out << part[list[i][j]].Position() << '\n';
   }
   else if(tp=="Direct")
   {
    Vector tmp = part[list[i][j]].Position();
    double lx = cell[0].Module();
    double ly = cell[1].Module();
    double lz = cell[2].Module();
    tmp = Vector(tmp[0]/lx,tmp[1]/ly,tmp[2]/lz);
    out << "  " << tmp << '\n';
   }
  }
 }
 
 // Print velocities, if there are.
 if (level==1)
 {
  out << '\n';
  for(long i=0 ; i < list.Size() ; ++i)
  {
   for (long int j=0;j<list[i].Size();++j)
   {
    out << "  " << part[list[i][j]].Velocity() << '\n';
   }
  }
 }

}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new VaspFormat(args); }
void destroy(Plugin * m) { delete m; }

