//
//
//


#include <sstream>

#include <lpmd/matrix.h>
#include <lpmd/util.h>

#include "angdist.h"

using namespace lpmd;

AngDist::AngDist(std::string args): Module("angdist")
{
 m = NULL;
 do_average = false;
 ProcessArguments(args);
}

AngDist::~AngDist()
{
 if (m != NULL) delete m;
}

//
void AngDist::SetParameter(std::string name)
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

void AngDist::Show() const
{
 Module::Show();
 std::cout << "   Atom Number = " << na << '\n';
 std::cout << "   Bins        = " << nb << '\n';
 std::cout << "   Atoms = ";
 for(unsigned int i=0;i<satoms.size();i++) std::cout << satoms[i] << "\t";
 std::cout << '\n';
 std::cout << "   Cutoffs = " << std::endl;
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
 std::cout << "   start = " << start_step << '\n';
 std::cout << "   end = " << end_step << '\n';
 std::cout << "   step = " << interval << '\n';
 std::cout << "   output = " << outputfile << '\n';
 std::cout << "   average = " << std::boolalpha << do_average << '\n';
}

std::string AngDist::Keywords() const { return "atoms rcut bins start end each output average"; }


void AngDist::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb <= 0 ||  na <=0)
 {
  std::cerr << "Error with ANG Calculation!!!." << std::endl;
 }

 int nsp = na;
 int N = simcell.Size();
 double **ang;
 ang = new double*[nb];
 for(int i=0;i<nb;i++) { ang[i]=new double[(int)(nsp*nsp*nsp)]; }
 for (int i=0;i<nb;i++) 
 { 
  for (int j=0;j<(int)(nsp*nsp*nsp);j++) ang[i][j]=0.0e0;
 }
 const std::list<std::string> lst = simcell.SpeciesTriplets();

 int s=0;
 for(std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  //Hace funcional el triplete de atomos y asigna los dos rcut
  std::vector<std::string> loa = SplitTextLine(*it,'-');
  int e1 = ElemNum(loa[0]);
  int e2 = ElemNum(loa[1]);
  int e3 = ElemNum(loa[2]);
  double rc12 = rcut[loa[0]+"-"+loa[1]];
  double rc23 = rcut[loa[1]+"-"+loa[2]];
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
       for(int k=0;k<N;++k)
       {
	if(simcell.GetAtom(k).Species()==e3)
	{
	 double distance2=simcell.Distance2(k,j);
	 if((k!=i && k!=j) && distance2<=rc23*rc23)
	 {
	  double angle=simcell.Angle(i,j,k);
	  int ig=(long)floor(nb*(angle/M_PI));
	  if(simcell.GetAtom(i).Species()==simcell.GetAtom(k).Species())
	  {
	   ang[ig][s]=ang[ig][s]+0.5;
  	  }
  	  else if(simcell.GetAtom(i).Species()!=simcell.GetAtom(k).Species())
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
  }
  s++; 
 }

 //
 // Output of ang
 //
 if (m != NULL) delete m;
 m = new Matrix(1 + nsp*nsp*nsp, nb);

 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "angle");
 int j=1;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, *it);
  j++;
 }
 //
 for(int i=0;i<nb;i++)
 {
  m->Set(0, i, (double)180*i/nb);
  for(int j=0;j<(int)(nsp*nsp*nsp);j++)
  {
   m->Set(j+1, i, ang[i][j]);
  }
 }
 for (int i=0;i<nb;i++) delete [] ang[i];
 delete [] ang;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new AngDist(args); }
void destroy(Module * m) { delete m; }

