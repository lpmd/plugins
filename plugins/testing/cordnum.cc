//
//
//

#include "cordnum.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/simulation.h>

#include <sstream>

using namespace lpmd;

CordNum::CordNum(std::string args): Plugin("cordnum", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("atoms");
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("rcut");
 DefineKeyword("maxn");
 DefineKeyword("output");
 DefineKeyword("average", "false");
 DefineKeyword("cutoff","0");
 ProcessArguments(args);
 nb = int(params["maxn"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 do_average = bool(params["average"]);
 cutoff = double(params["cutoff"]);
}

CordNum::~CordNum()
{
}

void CordNum::SetParameter(std::string name)
{
 #warning "SetParameter es horrible! mejorar los parametros"
 if (name == "atoms")
 {
  AssignParameter("atoms", GetNextWord());
  na = int((*this)["atoms"]);
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

void CordNum::Show(std::ostream & os) const
{
 Module::Show(os);
 std::cout << "   Atoms N     = " << na << '\n';
 std::cout << "   Atoms       = ";
 for(unsigned int i=0;i<satoms.size();i++) std::cout << satoms[i] << "\t";
 std::cout << std::endl;
 std::cout << "   Max neigh   = " << nb << '\n';
 std::cout << "   Cutoffs     = " << '\n';
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
   std::cout <<"\t"<< spec1 << "-" << spec2 << " = " << cut << std::endl;
  }
 }
}

void CordNum::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para calcular el numero de cordinacion de una     \n";
 std::cout << " celda de simulacion, utilizando los radios de corte entregados por el usuario.\n";
 std::cout << "      Se calcula el numero de cordinacion entre los vecinos de las subceldas   \n";
 std::cout << " generadas con el metodo linkedcell de la API liblpmd.                         \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      maxn          : Numero maximo de vecinos para escribir un histograma,    \n";
 std::cout << "                      usualmente un valor mayor 10.                            \n";
 std::cout << "      atoms         : Especifica el Numero de especies atomicas y sus simbolos \n";
 std::cout << "                      para el calculo de el numero de cordinacion              \n";
 std::cout << "      rcut          : Se especifican dos especies atomicas seguidas por su     \n";
 std::cout << "                      radio de corte.                                          \n";
 std::cout << "      output        : Archivo de salida para la informacion de la distribucion.\n";
 std::cout << "      average       : True/False Para promediar o no las distribuciones.       \n";
 std::cout << "      cutoff        : Radio de corte para la lista de vecinos, si no se asigna \n";
 std::cout << "                      un valor, toma el valor de 2*rcut                        \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cordnum                                                                   \n";
 std::cout << "     maxn 15                                                                   \n";
 std::cout << "     atoms 1 Ar                                                                \n";
 std::cout << "     rcut Ar Ar 3.95                                                           \n";
 std::cout << "     output cordnum.dat                                                        \n";
 std::cout << "     average false                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " property cordnum start=0 each=1 end=100                                     \n\n";
 std::cout << "      De esta forma calculamos el numeo de cordinacion de nuestra celda cada un\n";
 std::cout << " paso entre los pasos 0 y 100 de la simulacion de lpmd.                        \n";
}

void CordNum::Evaluate(Configuration & con, Potential & pot)
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::Array <int> esp = atoms.Elements();
 if (nb <= 0 || na <=0) throw PluginError("cordnum", "Error in coordination number calculation.");
 int nsp = na;
 unsigned long int N = atoms.Size();
 int **histo;
 double **cnfun;
 histo = new int*[nsp*nsp];
 cnfun = new double*[nsp*nsp];
 for(int i=0;i<nsp*nsp;i++) { histo[i]=new int[N]; }
 for(int i=0;i<nsp*nsp;i++) { cnfun[i]=new double[nb];}
 for (int i=0;i<nsp*nsp;i++)
 { 
  for (unsigned long int j=0;j<N;j++) histo[i][j]=0;
  for (int j=0;j<nb;j++) cnfun[i][j]=0.0e0;
 }

 //const std::list<std::string> lst = simcell.RepeatedSpeciesPairs(); 
 lpmd::Array<std::string> pairs;
 for(int i=0; i<esp.Size() ; ++i)
 {
  for(int j=0; j<esp.Size() ; ++j)
  {
   std::ostringstream ostr;
   ostr << ElemSym[esp[i]]<< "-" << ElemSym[esp[j]];
   pairs.Append(ostr.str());
  }
 }
 int s=0;

 for(int i=0;i<pairs.Size();++i)
 {
  //Hace funcional cada una de las especies de los pares.
  lpmd::Array<std::string> loa = StringSplit(pairs[i],'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  double rc12 = rcut[loa[0]+"-"+loa[1]];
  //Cuenta los atomos de la especie 1.
  int ne1=0;
  for (unsigned long int i=0;i<N;++i) {if(atoms[i].Z()==e1) ne1++;}
  //Comienzan las iteraciones.
  if (cutoff == 0) cutoff = rc12*2;
  for (unsigned long int i=0;i<N;++i)
  {
   if(atoms[i].Z()==e1)
   {
    lpmd::NeighborList & nlist = con.Neighbors(i,true,cutoff);
    for(long int k=0;k<nlist.Size();++k)
    {
     const lpmd::AtomPair & nn = nlist[k];
     if(nn.j->Z()==e2)
     {
      if(nn.r<=rc12)
      {
       histo[s][i]++;
      }
     }
    }
   }
  }
  for (unsigned long int i=0;i<N;i++)
  {
   if(atoms[i].Z()==e1 && histo[s][i]<nb)
   {
    cnfun[s][histo[s][i]]++;
   }
   histo[s][i]=0;
  }
  double distp=0.0e0;
  for(int i=0;i<nb;i++)
  {
   distp+=i*cnfun[s][i];
   cnfun[s][i]=(cnfun[s][i]/ne1);
  }
  s++;
 }
 //
 // Output of cordnum - histogram format
 //
 lpmd::Matrix & m = CurrentValue();
 m = lpmd::Matrix(1 + nsp*nsp, nb);
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "numb of neigh");
 int j=1;
 for (int i=0;i<pairs.Size();++i)
 {
  m.SetLabel(j, pairs[i]+" cn");
  j++;
 }
 //
 for(int i=0;i<nb;i++)
 {
  m.Set(0, i, i);
  for(int j=0;j<(int)(nsp*nsp);j++)
  {
   m.Set(j+1, i, cnfun[j][i]);
  }
 }
 //Borra arreglos dinamicos.
 for(int i=0;i<nsp*nsp;i++)
 {
  delete[] cnfun[i];
  delete[] histo[i];
 }
 delete [] cnfun;
 delete [] histo;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new CordNum(args); }
void destroy(Plugin * m) { delete m; }

