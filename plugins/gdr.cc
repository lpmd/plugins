//
//
//


#include <sstream>

#include <lpmd/matrix.h>
#include <lpmd/util.h>

#include "gdr.h"

using namespace lpmd;

Gdr::Gdr(std::string args): Module("gdr")
{
 m = NULL;
 do_average = false;
 ProcessArguments(args);
}

Gdr::~Gdr()
{
 if (m != NULL) delete m;
}

void Gdr::SetParameter(std::string name)
{
 if (name == "rcut") 
 {
  AssignParameter("rcut", GetNextWord());
  rcut = GetDouble("rcut");
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

void Gdr::Show() const
{
 Module::Show();
 std::cout << "   rcut     = " << rcut << '\n';
 std::cout << "   bins     = " << nb << '\n';
 std::cout << "   start    = " << start_step << '\n';
 std::cout << "   end      = " << end_step << '\n';
 std::cout << "   step     = " << interval << '\n';
 std::cout << "   output   = " << outputfile << '\n';
 std::cout << "   average  = " << std::boolalpha << do_average << '\n';
}

void Gdr::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = GDR " << '\n';
 std::cout << " Module Version     = 1.0 " << '\n';
 std::cout << " Support API lpmd   = 1.0 " << '\n';
 std::cout << " Problems Report to = gnm@gnm.cl " << '\n';
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Information "<< '\n' << '\n';
 std::cout << "      Modulo utilizado para calcular la funcion de autocorrelacion de pares \n";
 std::cout << " utiliza las condiciones de borde periodicas de la celda. \n";
 std::cout << "  " << '\n';
 std::cout << "      Los niveles de XYZ son 0,1 y 2 donde 0 muestra posiciones \n";
 std::cout << " 1 muestra posiciones y velocidades, 2 incluye aceleraciones\n\n";
 std::cout << " Use to read    = input module=xyz file=file-input.xyz level=0/1/2 \n" ;
 std::cout << " Use to write   = output module=xyz file=file-output.xyz each=steps level=0/1/2\n" ;
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example " << '\n' << '\n';
 std::cout << " Escribiendo    = output module=xyz file=fileoutput.xyz each=10 level=0\n" << '\n';
 std::cout << " En este caso el fichero de salida es de la forma xyz y se graba\n";
 std::cout << " cada 10 steps con un nivel 0" <<'\n';
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string Gdr::Keywords() const { return "rcut bins start end step output average"; }

void Gdr::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb == 0 || fabs(rcut) < 1e-05)     // fabs(rcut) < 1e-05 used to avoid comparing doubles
 {
  std::cerr << "Error with GDR Calculation!!!. Cutoff and Bins with wrong value." << std::endl;
 }
 double dr = rcut/ double(nb);

 int nsp = simcell.NEspec();
 int N = simcell.Size();
 double **g, *gt;
 g = new double*[nb];
 for(int i=0;i<nb;i++) { g[i]=new double[(int)(nsp*(nsp+1)/2)]; }
 gt = new double[nb]; //total gdr
 for (int i=0;i<nb;i++) 
 { 
  gt[i]=0.0e0;
  for (int j=0;j<(int)(nsp*(nsp+1)/2);j++) g[i][j]=0.0e0;
 }
 int s=0;
 const std::list<std::string> lst = simcell.SpeciesPairs();

 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)	   
 {
  //Hace funcional cada una de las especies de los pares.
  std::vector<std::string> loa = SplitTextLine(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  //Cuenta los atomos de cada especie atomica.
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(simcell.GetAtom(m).Species()==e1) ne1++;
   if(simcell.GetAtom(m).Species()==e2) ne2++;
  }
  //Comienza la iteracion principal para el calculo de g(r).
  double fac=1.0e0;
  if(e1==e2) fac=2.0e0;
  for(int i=0;i<N-1;i++)
  {
   if(simcell.GetAtom(i).Species()==e1)
   {
    for(int j=i+1;j<N;j++)
    {
     if(simcell.GetAtom(j).Species()==e2)
     {
      double distance=simcell.Distance2(i,j);
      if(distance<=rcut*rcut)
      {
       int ig=(long)floor(simcell.Distance(i,j)/dr);
       g[ig][s]=g[ig][s]+fac*(simcell.Volume())/(4.0e0*M_PI*simcell.Distance2(i,j)*dr*ne1*ne2);
      }
     }
    }
   }
  }
  s++;
 }
 //Calcula el valor de g(r) total.
 int j=0;
 for(std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  //Hace funcional cada una de las especies de los pares.
  std::vector<std::string> loa = SplitTextLine(*it,'-'); // lista de atomos
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  //Cuenta la concentracion atomica de cada especie atomica.
  int ne1=0,ne2=0;
  for(int m=0;m<N;m++)
  {
   if(simcell.GetAtom(m).Species()==e1) ne1++;
   if(simcell.GetAtom(m).Species()==e2) ne2++;
  }
  double ce1 = (double)ne1/(double)N;
  double ce2 = (double)ne2/(double)N;
  //Comienza la asignacion principal para g(r) total.
  for(int i=0;i<nb;i++)
  {
   if(e1==e2) gt[i] = gt[i]+ce1*ce2*g[i][j];
   else {gt[i]=gt[i]+2*ce1*ce2*g[i][j];}
  }
  j++;
 }
 //
 // Output of g(r)
 //
 if (m != NULL) delete m;
 m = new Matrix(2 + nsp*(nsp+1)/2, nb);

 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "r");
 m->SetLabel(nsp*(nsp+1)/2+1, "total g(r)");
 j=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, (*it)+" g(r)");
  j++;
 }
 // 
 for(int i=0;i<nb;i++)
 {
  m->Set(0, i, dr*i);
  for(int j=0;j<(int)(nsp*(nsp+1)/2);j++)
  {
   m->Set(j+1, i, g[i][j]);
  }
  m->Set(nsp*(nsp+1)/2+1, i, gt[i]);
 }
 delete [] gt;
 for (int i=0;i<nb;i++) delete [] g[i];
 delete [] g;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new Gdr(args); }
void destroy(Module * m) { delete m; }


