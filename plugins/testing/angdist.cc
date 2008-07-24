//
//
//


#include "angdist.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>

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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = angdist                                                  \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string AngDist::Keywords() const { return "atoms rcut bins start end each output average"; }

void AngDist::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (nb <= 0 ||  na <=0) throw PluginError("angdist", "Error with angular distribution calculation.");
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
   if(simcell[i].Species()==e2)
   {
    std::list<Neighbor> nlist;
    simcell.BuildNeighborList(i,nlist,true, simcell.CMCutoff());
    for(std::list<Neighbor>::const_iterator it=nlist.begin() ; it!=nlist.end() ; ++it)
    {
     const Neighbor & nn = *it;
     if(nn.j->Species()==e1 && nn.r*nn.r <= rc12*rc12)
     {
      for(std::list<Neighbor>::const_iterator jt=nlist.begin() ; jt!=nlist.end() ; ++jt)
      {
       const Neighbor & mm = *jt;
       if(&mm!=&nn)
       {
	if(mm.j->Species()==e3 && mm.r*mm.r <= rc23*rc23)
	{
	 double angle=simcell.Angle(nn.j->Index(),i,mm.j->Index());
	 int ig=(long)floor(nb*(angle/M_PI));
	 if(nn.j->Species()==mm.j->Species())
	 {
	  ang[ig][s]=ang[ig][s]+0.5;
	 }
	 else if(nn.j->Species()!=mm.j->Species())
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
