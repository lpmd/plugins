//
//
//

#include "cordnumfunc.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/simulationcell.h>

#include <sstream>

using namespace lpmd;

CordNumFunc::CordNumFunc(std::string args): Module("cordnumfunc")
{
 m = NULL;
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("atoms");
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("rcut");
 DefineKeyword("output");
 DefineKeyword("bins", "200");
 DefineKeyword("average", "false");
 ProcessArguments(args);
 cut = GetDouble("rcut");
 nb = GetInteger("bins");
 start = GetInteger("start");
 end = GetInteger("end");
 each = GetInteger("each");
 outputfile = GetString("output");
 do_average = GetBool("average");
}

CordNumFunc::~CordNumFunc()
{
 if (m != NULL) delete m;
}

void CordNumFunc::SetParameter(std::string name)
{
 if (name == "atoms")
 {
  AssignParameter("atoms", GetNextWord());
  na = GetInteger("atoms");
  for(int i=0;i<na;i++) { satoms.push_back(GetNextWord()); }
 }
 else Module::SetParameter(name);
}

void CordNumFunc::Show(std::ostream & os) const
{
 Module::Show(os);
 os << "   Atoms Number= " << na << '\n';
 os << "   Atoms       = ";
 for(unsigned int i=0;i<satoms.size();i++) std::cout << satoms[i] << "\t";
 os << std::endl;
}

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
 std::cout << "      atoms         : Especifica el Numero de especies atomicas y sus simbolos \n";
 std::cout << "                      para el calculo de el numero de cordinacion              \n";
 std::cout << "      rcut          : Se especifica el valor del radio de corte maximo de la   \n";
 std::cout << "                      funcion.                                                 \n";
 std::cout << "      output        : Archivo de salida para la informacion de la distribucion.\n";
 std::cout << "      average       : True/False Para promediar o no las distribuciones.       \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cordnumfunc                                                               \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     atoms 2 Ge O                                                              \n";
 std::cout << "     rcut 10.0                                                                 \n";
 std::cout << "     output cordnumfunc.dat                                                    \n";
 std::cout << "     average false                                                             \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " property cordnumfunc start=0 each=1 end=100                                 \n\n";
 std::cout << "      De esta forma calculamos el numero de cordinacion de nuestra celda cada  \n";
 std::cout << " un paso entre los pasos 0 y 100 de la simulacion de lpmd.                     \n";
}

void CordNumFunc::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb <= 0 || na <=0) throw PluginError("cordnumfunc", "Error in calculation.");
 int nsp = na;
 unsigned long int N = simcell.size();
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

 const std::list<std::string> lst = simcell.RepeatedSpeciesPairs();
 int s=0;
 for(std::list<std::string>::const_iterator it = lst.begin();it!=lst.end();++it)	   
 {
  //Hace funcional cada una de las especies de los pares.
  std::vector<std::string> loa = StringSplit< std::vector<std::string> >(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]); 
  //Cuenta los atomos de la especie 1.
  int ne1=0;
  for (unsigned long int i=0;i<N;i++) {if (simcell[i].Species()==e1) ne1++;}
  //Comienzan las iteraciones.
  for (unsigned long int i=0;i<N;i++)
  {
   if(simcell[i].Species()==e1)
   {
    std::vector<Neighbor> nlist;
    simcell.BuildNeighborList(i,nlist,true, rcut);
    for(unsigned long int k=0;k<nlist.size();++k)
    {
     const Neighbor &nn = nlist[k];
     if(nn.j->Species()==e2)
     {
      if(nn.r*nn.r<rcut*rcut)
      {
       int in=(long)floor(nn.r/dr);
       histo[s][in]++;
      }
     }
    }
   }
  }
  for(int i=0;i<nb;++i)
  {
   for(int j=0;j<i;j++)
   {
    cnfun[s][i]+=(histo[s][j]/ne1);
   }
  }
  s++;
 }

 //
 // Output of cordnum - histogram format
 //
 m = new Matrix(1 + nsp*nsp, nb);
 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "r");
 int j=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, (*it)+" cn");
  j++;
 }
 //
 double k=0.0e0;
 for(int i=0;i<nb;i++)
 {
  m->Set(0, i, k);
  for(int j=0;j<(int)(nsp*nsp);j++)
  {
   m->Set(j+1, i, cnfun[j][i]);
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
Module * create(std::string args) { return new CordNumFunc(args); }
void destroy(Module * m) { delete m; }

