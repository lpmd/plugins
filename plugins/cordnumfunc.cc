//
//
//


#include <sstream>

#include <lpmd/matrix.h>
#include <lpmd/util.h>

#include "cordnumfunc.h"

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

void CordNumFunc::Show() const
{
 Module::Show();
 std::cout << "   Atoms Number = " << na << '\n';
 std::cout << "   Atoms        = ";
 for(unsigned int i=0;i<satoms.size();i++) std::cout << satoms[i] << "\t";
 std::cout << std::endl;
 std::cout << "   Bins Number  = " << nb << '\n';
 std::cout << "   Cutoff       = " << cut << '\n';
 std::cout << "   start        = " << start_step << '\n';
 std::cout << "   end          = " << end_step << '\n';
 std::cout << "   step         = " << interval << '\n';
 std::cout << "   output       = " << outputfile << '\n';
 std::cout << "   average      = " << std::boolalpha << do_average << '\n';
}

void CordNumFunc::ShowHelp() const
{
}

std::string CordNumFunc::Keywords() const { return "atoms rcut bins start end each output average"; } 

void CordNumFunc::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb <= 0 ||  na <=0)
 {
  std::cerr << "Error with Coordination Number Calculation!!!." << std::endl;
 }

 int nsp = na;
 int N = simcell.Size();
 int **histo;
 double **cnfun;
 histo = new int*[nsp*nsp];
 cnfun = new double*[nsp*nsp];
 for(int i=0;i<nsp*nsp;i++) { histo[i]=new int[N]; }
 for(int i=0;i<nsp*nsp;i++) { cnfun[i]=new double[nb];}
 for (int i=0;i<nsp*nsp;i++)
 { 
  for (int j=0;j<N;j++) histo[i][j]=0;
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
  int dummy=0;
  for(double k=0.0;k<rcut;k+=dr)
  {
   for(int i=0;i<N;i++)
   {
    if(simcell.GetAtom(i).Species()==e1)
    {
     for(int j=0;j<N;j++)
     {
      if(j!=i && simcell.GetAtom(j).Species()==e2)
      {
       double distance=simcell.Distance2(i,j);
       if(distance >= k*k && distance < (k+dr)*(k+dr))
       {
	histo[s][i]++;
       }
      }
     }
    }
   }
   for(int i=0;i<N;i++)
   {
    if(simcell.GetAtom(i).Species()==e1)
    {
     cnfun[s][dummy]+=histo[s][i];
    }
   }
   dummy++;
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

