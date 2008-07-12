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
 do_average = false;
 ProcessArguments(args);
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
  for(int i=0;i<na;i++) 
  {
     satoms.push_back(GetNextWord());
  }
 }
 if (name == "rcut") 
 {
  AssignParameter("rcut",GetNextWord());
  cut = GetDouble("rcut");
 }
 if (name == "bins")
 {
  AssignParameter("bins", GetNextWord());
  nb = GetInteger("bins");
 }
 if (name == "start")
 {
  AssignParameter("start", GetNextWord());
  start_step = GetInteger("start");
 }
 if (name == "end")
 {
  AssignParameter("end", GetNextWord());
  end_step = GetInteger("end");
 }
 if (name == "each")
 {
  AssignParameter("each", GetNextWord());
  interval = GetInteger("each");
 }
 if (name == "output")
 {
  AssignParameter("output", GetNextWord());
  outputfile = GetString("output");
 }
 if (name == "average")
 {
  AssignParameter("average", GetNextWord());
  do_average = GetBool("average");
 }
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = cordnumfunc                                              \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

}

std::string CordNumFunc::Keywords() const { return "atoms rcut bins start end each output average"; } 

void CordNumFunc::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb <= 0 || na <=0) throw PluginError("cordnumfunc", "Error in calculation.");
 int nsp = na;
 int N = simcell.Size();
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
  std::vector<std::string> loa = SplitTextLine(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]); 
  //Cuenta los atomos de la especie 1.
  int ne1=0;
  for(int i=0;i<N;i++) {if(simcell.GetAtom(i).Species()==e1) ne1++;}
  //Comienzan las iteraciones.
  for(int i=0;i<N;i++)
  {
   if(simcell[i].Species()==e1)
   {
    std::list<Neighbor> nlist;
    simcell.BuildNeighborList(i,nlist,true, rcut);
    for(std::list<Neighbor>::const_iterator it=nlist.begin();it!=nlist.end();++it)
    {
     const Neighbor &nn = *it;
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

