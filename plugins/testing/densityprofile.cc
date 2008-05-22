//
//
//

#include "densityprofile.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>

using namespace lpmd;

DensityProfile::DensityProfile(std::string args): Module("densityprofile")
{
 m = NULL;
 do_average = false;
 range[0][0]=0.0e0;range[0][1]=0.0e0;
 range[1][0]=0.0e0;range[1][1]=0.0e0;
 range[2][0]=0.0e0;range[2][1]=0.0e0;
 ProcessArguments(args);
}

DensityProfile::~DensityProfile()
{
 if (m != NULL) delete m;
}

void DensityProfile::SetParameter(std::string name)
{
 if (name == "axis") 
 {
  AssignParameter("axis", GetNextWord());
  std::string eje=GetString("axis");
  if((eje == "X" || eje == "x") || eje == "0") axis=0;
  else if((eje == "Y" || eje == "y") || eje == "1") axis=1;
  else if((eje == "Z" || eje == "z") || eje == "2") axis=2;
  else 
  {
   axis=-1;
   ShowWarning("plugin densityprofile", "Wrong setting of axis.");
  }
 }
 if (name == "bins")
 {
  AssignParameter("bins", GetNextWord());
  bins = GetInteger("bins");
 }
 if (name == "range")
 {
  int tmp=0;
  std::string eje=GetNextWord();
  if((eje == "X" || eje == "x") || eje == "0") tmp=0;
  else if((eje == "Y" || eje == "y") || eje == "1") tmp=1;
  else if((eje == "Z" || eje == "z") || eje == "2") tmp=2;
  else 
  {
   tmp=-1;
   ShowWarning("plugin densityprofile", "Wrong setting of axis range.");
  }
  std::string smin=GetNextWord();
  if(smin=="ALL" || smin == "all")
  {
   range[tmp][0]=0;
   range[tmp][1]=0;
  }
  else
  {
   std::string smax=GetNextWord();
   range[tmp][0]=atof(smin.c_str());
   range[tmp][1]=atof(smax.c_str());
  }
  if(range[tmp][0]>=range[tmp][1] &&  smin!="all" && smin!="ALL") ShowWarning("densityprofile", "min and max values are not consistent.");
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

void DensityProfile::Show(std::ostream & os) const
{
 Module::Show(os);
 os << "   axis         = " << axis << '\n';
 os << "   range       = X  [" <<range[0][0]<<","<<range[0][1]<<"]"<< '\n';
 os << "               = Y  [" <<range[1][0]<<","<<range[1][1]<<"]"<< '\n';
 os << "               = Z  [" <<range[2][0]<<","<<range[2][1]<<"]"<< '\n';
}

void DensityProfile::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = densityprofile                                           \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular el perfil de la densidad de la celda de   \n";
 std::cout << " simulacion, actualmente solo unidimensionalmente.                             \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      axis           : Especifica el eje en el que se realizara el calulo.      \n";
 std::cout << "      bins          : Especifica el numero de divisiones para el eje.          \n";
 std::cout << "      range         : Especifica el rango para calculo de densidad.            \n";
 std::cout << "      output        : Fichero en el que se graba la densidad.                  \n";
 std::cout << "      average       : Setea si calcula o no el promedio de cada evaluacion.    \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use densityprofile                                                            \n";
 std::cout << "     axis X                                                                     \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     range Y 10 20                                                             \n";
 std::cout << "     range Z all                                                               \n";
 std::cout << "     range X all                                                               \n";
 std::cout << "     output filedensity.out                                                    \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property densityprofile start=1 each=10 end=100                             \n\n";
 std::cout << "      De esta forma calculamos la funcion de distribucion radial de pares en   \n";
 std::cout << " la simulacion entre 1 y 100 cada 10 pasos.                                    \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string DensityProfile::Keywords() const { return "rcut bins start end step output average"; }

void DensityProfile::Evaluate(SimulationCell & simcell, Potential & pot)
{
 if (bins == 0) throw PluginError("densityprofile", "Error in calculation: Wrong value for \"bins\".");

 //Vectores base, celda de simulacion.
 if(range[0][0]==0 && range[0][1]==0) {range[0][0]=0;range[0][1]=(simcell.GetVector(0)).Mod();}
 if(range[1][0]==0 && range[1][1]==0) {range[1][0]=0;range[1][1]=(simcell.GetVector(1)).Mod();}
 if(range[2][0]==0 && range[2][1]==0) {range[2][0]=0;range[2][1]=(simcell.GetVector(2)).Mod();}

 if(range[0][0]==range[0][1]) throw PluginError("densityprofile", "Error in cell range in axis X.");
 if(range[1][0]==range[1][1]) throw PluginError("densityprofile", "Error in cell range in axis Y.");
 if(range[2][0]==range[2][1]) throw PluginError("densityprofile", "Error in cell range in axis Z.");

 Vector na = simcell.GetVector(0); na.Norm();
 Vector nb = simcell.GetVector(1); nb.Norm();
 Vector nc = simcell.GetVector(2); nc.Norm();

 Vector la = na*range[0][1]-na*range[0][0];
 Vector lb = nb*range[1][1]-nb*range[1][0];
 Vector lc = nc*range[2][1]-nc*range[2][0];

 double dr=0.0e0,vol=0.0e0;

 if(axis==0) {dr=la.Mod()/double(bins);}
 else if(axis==1) {dr=lb.Mod()/double(bins);}
 else if(axis==2) {dr=lc.Mod()/double(bins);}
 else {throw PluginError("densityprofile", "Error in axis setting to set 'dr'!.");}

 vol = fabs(Dot(lc,Crux(la,lb)));

 double dvol = vol/bins; //delta de volumen de cada rango.

 int nsp = simcell.NEspec();
 int N = simcell.Size();
 double **rho, *rhot;
 rho = new double*[bins];
 for(int i=0;i<bins;i++) { rho[i]=new double[(int)(nsp)]; }
 rhot = new double[bins]; //densidad total
 for (int i=0;i<bins;i++) 
 { 
  rhot[i]=0.0e0;
  for (int j=0;j<(int)(nsp);j++) rho[i][j]=0.0e0;
 }
 int s=0;
 const std::list<std::string> lst = simcell.SpeciesList();

 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)	   
 {
  //Asigna la especie correspondiente.
  int e = ElemNum(*it);
  //Cuenta los atomos de la especie e.
  int ne=0;
  for(int m=0;m<N;m++)
  {
   if(simcell.GetAtom(m).Species()==e) ne++;
  }
  //Comienza la iteracion principal para el calculo.
  for(int i=0;i<N;++i)
  {
   if(simcell[i].Species()==e)
   {
    //vemos la ubicacion atomica respecto a nuestra "rejilla".
    lpmd::Vector position = simcell[i].Position();
    double m = simcell[i].Mass();
    double x = position.GetX();
    double y = position.GetY();
    double z = position.GetZ();
    if(x>=range[0][0] && x<=range[0][1])
    {
     if(y>=range[1][0] && y<=range[1][1])
     {
      if(z>=range[2][0] && z<=range[2][1])
      {
       //Entonces esta dentro de nuestros "ranges"
       //ahora la ubicamos segun nuestro eje preferencial 'axis'.
       double pp=0.0e0;
       if(axis==0) pp = x;
       else if(axis==1) pp = y;
       else if(axis==2) pp = z;
       else ShowWarning("plugin densityprofile", "Bad calculation of densityprofile, check your 'axis' option.");
       int ir = (long) floor(pp/dr);
       rho[ir][s] += m/dvol;
      }
     }
    }
   }
  }
  s++;
 }
 //Calcula el valor de rho(r) total.
 int j=0;
 for(std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  //Comienza la asignacion principal para g(r) total.
  for(int i=0;i<bins;i++)
  {
   rhot[i] += rho[i][j];
  }
  j++;
 }
 //
 // Output of rho(r)
 //
 if (m != NULL) delete m;
 m = new Matrix(3 + nsp, bins);

 // Asigna los labels al objeto Matrix para cada columna
 m->SetLabel(0, "r");
 m->SetLabel(1, "t");
 m->SetLabel(nsp+2, "total rho(r)");
 j=2;
 for (std::list<std::string>::const_iterator it=lst.begin();it!=lst.end();++it)
 {
  m->SetLabel(j, (*it)+" rho(r)");
  j++;
 }
 // 
 for(int i=0;i<bins;i++)
 {
  m->Set(0, i, dr*i);
  m->Set(1, i, counter);
  for(int j=0;j<(int)(nsp);j++)
  {
   m->Set(j+2, i, rho[i][j]);
  }
  m->Set(nsp+2, i, rhot[i]);
 }
 delete [] rhot;
 for (int i=0;i<bins;i++) delete [] rho[i];
 delete [] rho;
 counter++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new DensityProfile(args); }
void destroy(Module * m) { delete m; }


