//
//
//

#include "cordnumfunc.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/simulation.h>

#include <sstream>

using namespace lpmd;

CordNumFunc::CordNumFunc(std::string args): Plugin("cordnumfunc", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("rcut","10");
 DefineKeyword("output","cordnumfunc-output.dat");
 DefineKeyword("bins", "200");
 DefineKeyword("average", "false");
 ProcessArguments(args);
 cut = double(params["rcut"]);
 nb = int(params["bins"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 do_average = bool(params["average"]);
}

CordNumFunc::~CordNumFunc() { }

void CordNumFunc::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para calcular el numero de cordinacion de una     \n";
 std::cout << " celda de simulacion, utilizando el metodo clasico, no el histograma.          \n";
 std::cout << "      Se calcula el numero de cordinacion entre los vecinos de las subceldas   \n";
 std::cout << " generadas con el metodo linkedcell de la API liblpmd.                         \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Numero de seccion en el que se divide el rango de la     \n";
 std::cout << "                      funcion entre 0 y rcut.                                  \n";
 std::cout << "      rcut          : Se especifica el valor del radio de corte maximo de la   \n";
 std::cout << "                      funcion.                                                 \n";
 std::cout << "      output        : Archivo de salida para la informacion de la distribucion.\n";
 std::cout << "      average       : True/False Para promediar o no las distribuciones.       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cordnumfunc                                                               \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     rcut 10.0                                                                 \n";
 std::cout << "     output cordnumfunc.dat                                                    \n";
 std::cout << "     average false                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " property cordnumfunc start=0 each=1 end=100                                 \n\n";
 std::cout << "      De esta forma calculamos el numero de cordinacion de nuestra celda cada  \n";
 std::cout << " un paso entre los pasos 0 y 100 de la simulacion de lpmd.                     \n";
}

void CordNumFunc::Evaluate(lpmd::Configuration & con, Potential & pot)
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::Array <int> esp = atoms.Elements();
 if (nb <= 0 || atoms.Size() <=0) throw PluginError("cordnumfunc", "Error in calculation.");
 int nsp = esp.Size();
 unsigned long int N = atoms.Size();
 double **histo;
 double **cnfun;
 histo = new double*[nsp*nsp];
 cnfun = new double*[nsp*nsp];
 for(int i=0;i<nsp*nsp;i++) { histo[i]=new double[nb]; }
 for(int i=0;i<nsp*nsp;i++) { cnfun[i]=new double[nb];}
 for (int i=0;i<nsp*nsp;i++)
 { 
  for (int j=0;j<nb;j++) histo[i][j]=0.0e0;
  for (int j=0;j<nb;j++) cnfun[i][j]=0.0e0;
 }
 double rcut=cut;
 double dr=rcut/nb;

 
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
  //Cuenta los atomos de la especie 1.
  int ne1=0;
  for (unsigned long int j=0;j<N;j++) {if (atoms[j].Z()==e1) ne1++;}
  //Comienzan las iteraciones.
  for (unsigned long int j=0;j<N;j++)
  {
   if(atoms[j].Z()==e1)
   {
    lpmd::NeighborList & nlist = con.Neighbors(j,true,rcut);
    for(long int k=0;k<nlist.Size();++k)
    {
     const lpmd::AtomPair & nn = nlist[k];
     if(nn.j->Z()==e2)
     {
      if(nn.r2<rcut*rcut)
      {
       int in=(long)floor(sqrt(nn.r2)/dr);
       histo[s][in]++;
      }
     }
    }
   }
  }
  for(int k=0;k<nb;++k)
  {
   for(int j=0;j<k;j++)
   {
    cnfun[s][k]+=(histo[s][j]/ne1);
   }
  }
  s++;
 }

 //
 // Output of cordnum - histogram format
 //
 Matrix & m = CurrentValue(); 
 m = lpmd::Matrix(1 + nsp*nsp, nb);
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "r");
 int j=1;
 for (int i=0;i<pairs.Size();++i)
 {
  m.SetLabel(j, pairs[i]+" cn");
  j++;
 }
 //
 double k=0.0e0;
 for(int i=0;i<nb;i++)
 {
  m.Set(0, i, k);
  for(j=0;j<(int)(nsp*nsp);j++)
  {
   m.Set(j+1, i, cnfun[j][i]);
  }
  k+=dr;
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
Plugin * create(std::string args) { return new CordNumFunc(args); }
void destroy(Plugin * m) { delete m; }

