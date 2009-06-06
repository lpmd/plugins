//
//
//


#include "angdist.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/atompair.h>
#include <lpmd/configuration.h>

#include <sstream>

using namespace lpmd;

AngDist::AngDist(std::string args): Plugin("angdist", "2.0")
{
 ParamList & param = (*this);
 m = NULL;
 //
 DefineKeyword("atoms");
 DefineKeyword("rcut");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("output");
 DefineKeyword("bins", "200");
 DefineKeyword("average", "false");
 // 
 ProcessArguments(args);
 nb = int(param["bins"]);
 start = int(param["start"]);
 end = int(param["end"]);
 each = int(param["each"]);
 OutputFile() = param["output"];
 do_average = bool(param["average"]);
}

AngDist::~AngDist()
{
 if (m != NULL) delete m;
}

void AngDist::SetParameter(std::string name)
{
 #warning "SetParameter es horrible! mejorar los parametros"
 if (name == "atoms")
 {
  AssignParameter("atoms", GetNextWord());
  na = (*this)["atoms"];
  for(int i=0;i<na;i++) { satoms.push_back(GetNextWord()); }
 }
 else if (name == "rcut") 
 {
  std::string atom1 = GetNextWord();
  std::string atom2 = GetNextWord();
  std::string tmp = GetNextWord();
  double cutoff = atof(tmp.c_str());
  rcut[atom1+"-"+atom2] = cutoff;
  rcut[atom2+"-"+atom1] = cutoff;
 }
 else Module::SetParameter(name);
}

void AngDist::Show(std::ostream & os) const
{
 Module::Show(os);
 os << "   Atom Number = " << na << '\n';
 os << "   Atoms       = ";
 for (unsigned int i=0;i<satoms.size();i++) os << satoms[i] << "\t";
 os << '\n';
 os << "   Cutoffs     = " << std::endl;
 for(unsigned int i=0;i<satoms.size();i++)
 {
  for(unsigned int j=i;j<satoms.size();j++)
  {
   std::string spec1=satoms[i];
   std::string spec2=satoms[j];
   std::string tmp = spec1+"-"+spec2;
   // Truco para acceder a rcut[tmp] desde un metodo const 
   const std::map<std::string, double>::const_iterator & p = rcut.find(tmp);
   double cut = (*p).second;
   // fin del truco
   os << "\t" << spec1 << "-" << spec2 << " = " << cut << std::endl;
  }
 }
}

void AngDist::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para calcular la distribucion angular de una celda\n";
 std::cout << " de simulacion, utilizando los radios de corte entregados por el usuario.      \n";
 std::cout << "      Se calcula la distribucion angular entre los vecinos de las subceldas    \n";
 std::cout << " generadas con el metodo linkedcell de la API liblpmd.                         \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Determina el numero de intervalos en los que se divide la\n";
 std::cout << "                      malla entre 0 y 180 grados.                              \n";
 std::cout << "      atoms         : Especifica el Numero de especies atomicas y sus simbolos \n";
 std::cout << "                      para el calculo de la distribucion angular               \n";
 std::cout << "      rcut          : Se especifican dos especies atomicas seguidas por su     \n";
 std::cout << "                      radio de corte.                                          \n";
 std::cout << "      output        : Archivo de salida para la informacion de la distribucion.\n";
 std::cout << "      average       : True/False Para promediar o no las distribuciones.       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use angdist                                                                   \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     atoms 2 Ge O                                                              \n";
 std::cout << "     rcut Ge Ge 3.60                                                           \n";
 std::cout << "     rcut Ge O  1.90                                                           \n";
 std::cout << "     rcut O  O  3.20                                                           \n";
 std::cout << "     output angdist.dat                                                        \n";
 std::cout << "     average false                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " property angdist start=0 each=1 end=100                                     \n\n";
 std::cout << "      De esta forma calculamos la distribucion angular de nuestra celda cada un\n";
 std::cout << " paso entre los pasos 0 y 100 de la simulacion de lpmd.                        \n";
}

void AngDist::Evaluate(lpmd::Configuration & con, lpmd::Potential & pot)
{
 lpmd::BasicParticleSet & part = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 if (nb <= 0 ||  na <=0) throw PluginError("angdist", "Error with angular distribution calculation.");
 int nsp = na;
 unsigned long int N = part.Size();
 double **ang;
 ang = new double*[nb];
 for(int i=0;i<nb;i++) { ang[i]=new double[(int)(nsp*nsp*nsp)]; }
 for (int i=0;i<nb;i++) 
 { 
  for (int j=0;j<(int)(nsp*nsp*nsp);j++) ang[i][j]=0.0e0;
 }
#warning Armado propio de triples, ¿debería ser método de la api?
 lpmd::Array <int> elements = part.Elements();
 Array<std::string> lst;
 for (int i=0;i<elements.Size();++i)
 {
  for (int j=0;j<elements.Size();++i)
  {
   for (int k=0;k<elements.Size();++i)
   {
    std::ostringstream ostr ;
    ostr << ElemSym[elements[i]] << "-" << ElemSym[elements[j]] << "-" <<ElemSym[elements[k]];
    lst.Append(ostr.str());
   }
  }
 }

 int s=0;
 for(int i=0;i<lst.Size();++i)
 {
  //Hace funcional el triplete de atomos y asigna los dos rcut
  Array<std::string> loa = StringSplit(lst[i],'-');
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  int e3 = ElemNum(loa[2]);
  double rc12 = rcut[loa[0]+"-"+loa[1]];
  double rc23 = rcut[loa[1]+"-"+loa[2]];
  for(unsigned long int i=0;i<N;++i)
  {
   if (part[i].Z()==e2)
   {
#warning Que pasará com CMCutoff() ??? ahora esta rempplazado por rc12+rc23
    lpmd::NeighborList & nlist = con.Neighbors(i,true,rc12+rc23);
    for (long int k=0;k<nlist.Size();++k)
    {
     const lpmd::AtomPair & nn = nlist[k];
     if(nn.j->Z()==e1 && nn.r*nn.r <= rc12*rc12)
     {
      for (long int l=0;l<nlist.Size();++l)
      {
       const lpmd::AtomPair & mm = nlist[l];
       if(&mm!=&nn)
       {
	if(mm.j->Z()==e3 && mm.r*mm.r <= rc23*rc23)
	{
#warning Angle debe pasarse a la API!
	 Vector a = cell.Displacement(part[i].Position(), nn.j->Position());
	 Vector b = cell.Displacement(part[i].Position(), mm.j->Position());
	 double angle= acos(Dot(a,b)/ (a.Module()*b.Module()));
	 int ig=(long)floor(nb*(angle/M_PI));
	 if(nn.j->Z()==mm.j->Z())
	 {
	  ang[ig][s]=ang[ig][s]+0.5;
	 }
	 else if(nn.j->Z()!=mm.j->Z())
	 {
	  ang[ig][s]=ang[ig][s]+1.0;
	 }
	}
       }
      }
     }
    }
   }
  }
  s++; 
 }

 //
 // Output of ang
 //
 lpmd::Matrix & m = CurrentValue();
 m = lpmd::Matrix(1 + nsp*nsp*nsp, nb);

 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "angle");
 int j=1;
 for (int i=0;i < lst.Size() ; ++i)
 {
  m.SetLabel(j, lst[i]);
  j++;
 }
 //
 for(int i=0;i<nb;i++)
 {
  m.Set(0, i, (double)180*i/nb);
  for(int j=0;j<(int)(nsp*nsp*nsp);j++)
  {
   m.Set(j+1, i, ang[i][j]);
  }
 }
 for (int i=0;i<nb;i++) delete [] ang[i];
 delete [] ang;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new AngDist(args); }
void destroy(Plugin * m) { delete m; }

