//
//
//

#include "tempprofile.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/simulation.h>
#include <lpmd/plugin.h>
#include <lpmd/session.h>

#include <sstream>

using namespace lpmd;

TempProfile::TempProfile(std::string args): Plugin("tempprofile", "2.0")
{
 ParamList & param = (*this);
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("rcut");
 DefineKeyword("output");
 DefineKeyword("bins", "200");
 DefineKeyword("average", "false");
 range[0][0]=0.0e0;range[0][1]=0.0e0;
 range[1][0]=0.0e0;range[1][1]=0.0e0;
 range[2][0]=0.0e0;range[2][1]=0.0e0;
 ProcessArguments(args);
 bins = int(param["bins"]);
 start = int(param["start"]);
 end = int(param["end"]);
 each = int(param["each"]);
 OutputFile() = param["output"];
 do_average = bool(param["average"]);
}

TempProfile::~TempProfile()
{
}

void TempProfile::SetParameter(std::string name)
{
 //#warning "SetParameter es horrible! mejorar los parametros"
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
   ShowWarning("plugin tempprofile", "Wrong setting of axis.");
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
   ShowWarning("plugin tempprofile", "Wrong setting of axis range.");
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
  if(range[tmp][0]>=range[tmp][1] &&  smin!="all" && smin!="ALL") ShowWarning("tempprofile", "min and max values are not consistent.");
 }
 else Module::SetParameter(name);
}

void TempProfile::Show(std::ostream & os) const
{
 Module::Show(os);
 os << "   axis         = " << axis << '\n';
 os << "   range       = X  [" <<range[0][0]<<","<<range[0][1]<<"]"<< '\n';
 os << "               = Y  [" <<range[1][0]<<","<<range[1][1]<<"]"<< '\n';
 os << "               = Z  [" <<range[2][0]<<","<<range[2][1]<<"]"<< '\n';
}

void TempProfile::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para calcular el perfil de temperaturas de la celda de  \n";
 std::cout << " simulacion, actualmente solo unidimensionalmente en el tiempo.                \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      axis          : Especifica el eje en el que se realizara el calulo.      \n";
 std::cout << "      bins          : Especifica el numero de divisiones para el eje.          \n";
 std::cout << "      range         : Especifica el rango para calculo de densidad.            \n";
 std::cout << "      output        : Fichero en el que se graba la densidad.                  \n";
 std::cout << "      average       : Setea si calcula o no el promedio de cada evaluacion.    \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use tempprofile                                                            \n";
 std::cout << "     axis X                                                                     \n";
 std::cout << "     bins 200                                                                  \n";
 std::cout << "     range Y 10 20                                                             \n";
 std::cout << "     range Z all                                                               \n";
 std::cout << "     range X all                                                               \n";
 std::cout << "     output filetemperature.out                                                \n";
 std::cout << "     average true                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property tempprofile start=1 each=10 end=100                                \n\n";
 std::cout << "      De esta forma calculamos la funcion de distribucion radial de pares en   \n";
 std::cout << " la simulacion entre 1 y 100 cada 10 pasos.                                    \n";
}

void TempProfile::Evaluate(Configuration & con, Potential & pot)
{
 lpmd::BasicParticleSet & atoms = con.Atoms();
 lpmd::BasicCell & cell = con.Cell();
 lpmd::Array<int> elements = atoms.Elements();
 if (bins == 0) throw PluginError("tempprofile", "Error in calculation: Wrong value for \"bins\".");

 //Vectores base, celda de simulacion.
 if(range[0][0]==0 && range[0][1]==0) {range[0][0]=0;range[0][1]=cell[0].Module();}
 if(range[1][0]==0 && range[1][1]==0) {range[1][0]=0;range[1][1]=cell[1].Module();}
 if(range[2][0]==0 && range[2][1]==0) {range[2][0]=0;range[2][1]=cell[2].Module();}

 if(range[0][0]==range[0][1]) throw PluginError("tempprofile", "Error in cell range in axis X.");
 if(range[1][0]==range[1][1]) throw PluginError("tempprofile", "Error in cell range in axis Y.");
 if(range[2][0]==range[2][1]) throw PluginError("tempprofile", "Error in cell range in axis Z.");

 Vector na = cell[0]; na.Normalize();
 Vector nb = cell[1]; nb.Normalize();
 Vector nc = cell[2]; nc.Normalize();

 Vector la = na*range[0][1]-na*range[0][0];
 Vector lb = nb*range[1][1]-nb*range[1][0];
 Vector lc = nc*range[2][1]-nc*range[2][0];

 double dr=0.0e0;

 if(axis==0) {dr=la.Module()/double(bins);}
 else if(axis==1) {dr=lb.Module()/double(bins);}
 else if(axis==2) {dr=lc.Module()/double(bins);}
 else {throw PluginError("tempprofile", "Error in axis setting to set 'dr'!.");}

 int nsp = elements.Size();
 long int N = atoms.Size();
 long int *nab = new long int[bins];
 double **temp, *tempt;
 temp = new double*[bins];
 for(int i=0;i<bins;i++) { temp[i]=new double[(int)(nsp)]; }
 tempt = new double[bins]; //temperatura total
 for (int i=0;i<bins;i++) 
 { 
  tempt[i]=0.0e0;
  nab[i] = 0;
  for (int j=0;j<(int)(nsp);j++) temp[i][j]=0.0e0;
 }
 int s=0;
 const double kin2ev = double(GlobalSession["kin2ev"]);
 const double kboltzmann = double(GlobalSession["kboltzmann"]);
 for (int i=0;i<elements.Size();++i)	   
 {
  //Asigna la especie correspondiente.
  int e = elements[i];
  //Cuenta los atomos de la especie e.
  int ne=0;
  for(long int m=0;m<N;m++)
  {
   if(atoms[m].Z()==e) ne++;
  }
  //Comienza la iteracion principal para el calculo.
  for(long int i=0;i<N;++i)
  {
   if(atoms[i].Z()==e)
   {
    //vemos la ubicacion atomica respecto a nuestra "rejilla".
    lpmd::Vector position = atoms[i].Position();
    lpmd::Vector velocity = atoms[i].Velocity();
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
       else ShowWarning("plugin tempprofile", "Bad calculation of tempprofile, check your 'axis' option.");
       int ir = (long) floor(pp/dr);
       temp[ir][s] += (1/2)*m*velocity.SquareModule()*kin2ev; //Solo son los aportes a la energia cinetica.
       nab[ir]++;
      }
     }
    }
   }
  }
  //Reasigna el verdadero valor de tempratura, segun la cantidad de atomos de la especie considerada 
  //en cada bins (nab).
  for(int i=0;i<bins;i++)
  {
   temp[i][s] = ((2.0/3.0)*temp[i][s])/(kboltzmann*double(nab[i]));
  }
  s++;
 }
 //Calcula el valor de temp(r) total.
 int j=0;
 for(int i=0;i<elements.Size();++i)
 {
  for(int i=0;i<bins;i++)
  {
   tempt[i] += temp[i][j];
  }
  j++;
 }
 //
 // Output of rho(r)
 //
 lpmd::Matrix & m  = CurrentValue();
 m = lpmd::Matrix(3 + nsp, bins);

 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "r");
 m.SetLabel(1, "t");
 m.SetLabel(nsp+2, "total Temp(r)");
 j=2;
 for (int i=0;i<elements.Size();++i)
 {
  m.SetLabel(j, elements[i]+" temp(r)");
  j++;
 }
 // 
 for(int i=0;i<bins;i++)
 {
  m.Set(0, i, dr*i);
  m.Set(1, i, counter);
  for(int j=0;j<(int)(nsp);j++)
  {
   m.Set(j+2, i, temp[i][j]);
  }
  m.Set(nsp+2, i, tempt[i]);
 }
 delete [] tempt;
 for (int i=0;i<bins;i++) delete [] temp[i];
 delete [] temp;
 delete [] nab; 
 counter++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new TempProfile(args); }
void destroy(Plugin * m) { delete m; }


