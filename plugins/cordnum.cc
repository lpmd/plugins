//
//
//


#include <sstream>

#include <lpmd/matrix.h>
#include <lpmd/util.h>

#include "cordnum.h"

using namespace lpmd;

CordNum::CordNum(std::string args): Module("cordnum")
{
 m = NULL;
 do_average = false;
 ProcessArguments(args);
}

CordNum::~CordNum()
{
 if (m != NULL) delete m;
}

void CordNum::SetParameter(std::string name)
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
  std::string atom1 = GetNextWord();
  std::string atom2 = GetNextWord();
  std::string tmp = GetNextWord();
  double cutoff = atof(tmp.c_str());
  rcut[atom1+"-"+atom2] = cutoff;
  rcut[atom2+"-"+atom1] = cutoff;
 }
 if (name == "maxn")
 {
  AssignParameter("maxn", GetNextWord());
  nb = GetInteger("maxn");
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

void CordNum::Show() const
{
 Module::Show();
 std::cout << "   Atoms N                  = " << na << '\n';
 std::cout << "   Atoms                    = ";
 for(unsigned int i=0;i<satoms.size();i++) std::cout << satoms[i] << "\t";
 std::cout << std::endl;
 std::cout << "   Max number of neighbours = " << nb << '\n';
 std::cout << "   Cutoffs                  = " << '\n';
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
 std::cout << "   start                    = " << start_step << '\n';
 std::cout << "   end                      = " << end_step << '\n';
 std::cout << "   step                     = " << interval << '\n';
 std::cout << "   output                   = " << outputfile << '\n';
 std::cout << "   average                  = " << std::boolalpha << do_average << '\n';
}

void CordNum::ShowHelp() const
{
}

std::string CordNum::Keywords() const { return "atoms rcut maxn start end each output average"; } 

void CordNum::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb <= 0 ||  na <=0){std::cerr << "Error with Coordination Number Calculation!!!." << std::endl;}

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

 const std::list<std::string> lst = simcell.RepeatedSpeciesPairs(); 
 int s=0;

 for(std::list<std::string>::const_iterator it = lst.begin() ; it!=lst.end() ; ++it)
 {
  //Hace funcional cada una de las especies de los pares.
  std::vector<std::string> loa = SplitTextLine(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]); 
  double rc12 = rcut[loa[0]+"-"+loa[1]];
  //Cuenta los atomos de la especie 1.
  int ne1=0;
  for(int i=0;i<N;i++) {if(simcell.GetAtom(i).Species()==e1) ne1++;}
  //Comienzan las iteraciones.
  for(int i=0;i<N;++i)
  {
   if(simcell.GetAtom(i).Species()==e1)
   {
    for(int j=0;j<N;++j)
    {
     if(j!=i && simcell.GetAtom(j).Species()==e2)
     {
      double distance=simcell.Distance2(i,j);
      if(distance<=rc12*rc12)
      {
       histo[s][i]++;
      }
     }
    }
   }	  
  }
  for(int i=0;i<N;i++)
  {
   if(simcell.GetAtom(i).Species()==e1 && histo[s][i]<nb)
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
 m = new Matrix(1 + nsp*nsp, nb);
 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "numb of neigh");
 int j=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, (*it)+" cn");
  j++;
 }
 //
 for(int i=0;i<nb;i++)
 {
  m->Set(0, i, i);
  for(int j=0;j<(int)(nsp*nsp);j++)
  {
   m->Set(j+1, i, cnfun[j][i]);
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
Module * create(std::string args) { return new CordNum(args); }
void destroy(Module * m) { delete m; }

