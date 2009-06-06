//
//
//

#include "densityprofile.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <lpmd/plugin.h>

#include <sstream>

using namespace lpmd;

DensityProfile::DensityProfile(std::string args): Plugin("densityprofile", "2.0")
{
 lpmd::ParamList & params = (*this);
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("rcut");
 DefineKeyword("output");
 DefineKeyword("bins", "200");
 DefineKeyword("average", "false");
 range[0][0]=0.0e0;range[0][1]=0.0e0;
 range[1][0]=0.0e0;range[1][1]=0.0e0;
 range[2][0]=0.0e0;range[2][1]=0.0e0;
 ProcessArguments(args);
 bins = int(params["bins"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 do_average = bool(params["average"]);
}

DensityProfile::~DensityProfile()
{
}

void DensityProfile::SetParameter(std::string name)
{
 #warning "SetParameter es horrible! mejorar los parametros"
 if (name == "axis") 
 {
  AssignParameter("axis", GetNextWord());
  std::string eje=(*this)["axis"];
  if((eje == "X" || eje == "x") || eje == "0") axis=0;
  else if((eje == "Y" || eje == "y") || eje == "1") axis=1;
  else if((eje == "Z" || eje == "z") || eje == "2") axis=2;
  else 
  {
   axis=-1;
   lpmd::ShowWarning("plugin densityprofile", "Wrong setting of axis.");
  }
 }
 else if (name == "range")
 {
  int tmp=0;
  std::string eje=GetNextWord();
  if((eje == "X" || eje == "x") || eje == "0") tmp=0;
  else if((eje == "Y" || eje == "y") || eje == "1") tmp=1;
  else if((eje == "Z" || eje == "z") || eje == "2") tmp=2;
  else 
  {
   tmp=-1;
   lpmd::ShowWarning("plugin densityprofile", "Wrong setting of axis range.");
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
  if(range[tmp][0]>=range[tmp][1] &&  smin!="all" && smin!="ALL") lpmd::ShowWarning("densityprofile", "min and max values are not consistent.");
 }
 else Module::SetParameter(name);
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
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular el perfil de la densidad de la celda de   \n";
 std::cout << " simulacion, actualmente solo unidimensionalmente.                             \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      axis           : Especifica el eje en el que se realizara el calulo.      \n";
 std::cout << "      bins          : Especifica el numero de divisiones para el eje.          \n";
 std::cout << "      range         : Especifica el rango para calculo de densidad.            \n";
 std::cout << "      output        : Fichero en el que se graba la densidad.                  \n";
 std::cout << "      average       : Setea si calcula o no el promedio de cada evaluacion.    \n";
 std::cout << '\n';
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
}

void DensityProfile::Evaluate(lpmd::Configuration & con, lpmd::Potential & pot)
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 lpmd::Array<int> elements = atoms.Elements();
 if (bins == 0) throw lpmd::PluginError("densityprofile", "Error in calculation: Wrong value for \"bins\".");

 //Vectores base, celda de simulacion.
 if(range[0][0]==0 && range[0][1]==0) {range[0][0]=0;range[0][1]=(cell[0]).Module();}
 if(range[1][0]==0 && range[1][1]==0) {range[1][0]=0;range[1][1]=(cell[1]).Module();}
 if(range[2][0]==0 && range[2][1]==0) {range[2][0]=0;range[2][1]=(cell[2]).Module();}

 if(range[0][0]==range[0][1]) throw lpmd::PluginError("densityprofile", "Error in cell range in axis X.");
 if(range[1][0]==range[1][1]) throw lpmd::PluginError("densityprofile", "Error in cell range in axis Y.");
 if(range[2][0]==range[2][1]) throw lpmd::PluginError("densityprofile", "Error in cell range in axis Z.");

 lpmd::Vector na = cell[0]; na.Normalize();
 lpmd::Vector nb = cell[1]; nb.Normalize();
 lpmd::Vector nc = cell[2]; nc.Normalize();

 lpmd::Vector la = na*range[0][1]-na*range[0][0];
 lpmd::Vector lb = nb*range[1][1]-nb*range[1][0];
 lpmd::Vector lc = nc*range[2][1]-nc*range[2][0];

 double dr=0.0e0,vol=0.0e0;

 if(axis==0) {dr=la.Module()/double(bins);}
 else if(axis==1) {dr=lb.Module()/double(bins);}
 else if(axis==2) {dr=lc.Module()/double(bins);}
 else {throw lpmd::PluginError("densityprofile", "Error in axis setting to set 'dr'!.");}

 vol = fabs(Dot(lc,Cross(la,lb)));

 double dvol = vol/bins; //delta de volumen de cada rango.

 int nsp = elements.Size();
 unsigned long int N = atoms.Size();
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

 for (int i=0;i<elements.Size();++i)	   
 {
  //Asigna la especie correspondiente.
  //Cuenta los atomos de la especie e.
  int ne=0;
  for(unsigned long int m=0;m<N;m++)
  {
   if(atoms[m].Z()==elements[i]) ne++;
  }
  //Comienza la iteracion principal para el calculo.
  for(unsigned long int i=0;i<N;++i)
  {
   if(atoms[i].Z()==elements[i])
   {
    //vemos la ubicacion atomica respecto a nuestra "rejilla".
    lpmd::Vector position = atoms[i].Position();
    double m = atoms[i].Mass();
    double x = position[0];
    double y = position[1];
    double z = position[2];
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
       else lpmd::ShowWarning("plugin densityprofile", "Bad calculation of densityprofile, check your 'axis' option.");
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
 for(int i=0;i<elements.Size();++i)
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
 lpmd::Matrix & m = CurrentValue();
 m = lpmd::Matrix(3 + nsp, bins);

 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "r");
 m.SetLabel(1, "t");
 m.SetLabel(nsp+2, "total rho(r)");
 j=2;
 for (int i=0;i<elements.Size();++i)
 {
  m.SetLabel(j, elements[i]+" rho(r)");
  j++;
 }
 // 
 for(int i=0;i<bins;i++)
 {
  m.Set(0, i, dr*i);
  m.Set(1, i, counter);
  for(int j=0;j<(int)(nsp);j++)
  {
   m.Set(j+2, i, rho[i][j]);
  }
  m.Set(nsp+2, i, rhot[i]);
 }
 delete [] rhot;
 for (int i=0;i<bins;i++) delete [] rho[i];
 delete [] rho;
 counter++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new DensityProfile(args); }
void destroy(Plugin * m) { delete m; }


