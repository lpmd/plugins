//
//
//

#include "veldist.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/configuration.h>

#include <sstream>

using namespace lpmd;

VelDist::VelDist(std::string args): Module("veldist")
{
 ParamList & params = (*this);
 AssignParameter("version", "1.0"); 
 AssignParameter("apirequired", "1.1"); 
 AssignParameter("bugreport", "gnm@gnm.cl"); 
 //
 DefineKeyword("start");
 DefineKeyword("end");
 DefineKeyword("each");
 DefineKeyword("output");
 DefineKeyword("bins", "300");
 ProcessArguments(args);
 bins = int(params["bins"]);
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
}

VelDist::~VelDist() { }

void VelDist::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "   Calcula el histograma de velocidades, tanto por magnitud como por componentes.\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      bins          : Especifica el numero de intervalos (bins) del histograma \n";
 std::cout << "      output        : Fichero en el que se graba la salida                     \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use veldist                                                                   \n";
 std::cout << "     output veldist.dat                                                        \n";
 std::cout << "     bins 500                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property veldist start=1 each=10 end=100                                      \n\n";
}

void VelDist::Evaluate(Configuration & conf, Potential & pot)
{
 //
 // Primero determina los valores extremos
 // 
 BasicParticleSet & atoms = conf.Atoms();
 long int n = atoms.Size();
 Vector vmin, vmax;
 double vmodmin, vmodmax, vrmin, vrmax;
 vmin = vmax = atoms[0].Velocity();
 vmodmin = vmodmax = atoms[0].Velocity().Module();
 vrmin = vmodmin;
 vrmax = vmodmax;
 for (long int i=0;i<n;++i)
 {
  const Vector & v = atoms[i].Velocity();
  for (int q=0;q<3;++q)
  {
   if (v[q] < vmin[q]) vmin[q] = v[q]; 
   if (v[q] > vmax[q]) vmax[q] = v[q]; 
   if (vmin[q] < vrmin) vrmin = vmin[q];
   if (vmax[q] > vrmax) vrmax = vmax[q];
  }
  if (v.Module() < vmodmin) vmodmin = v.Module();
  if (v.Module() > vmodmax) vmodmax = v.Module();
  if (vmodmin < vrmin) vrmin = vmodmin;
  if (vmodmax > vrmax) vrmax = vmodmax;
 }
 // 
 // Ahora llena el histograma 
 //
 Matrix & m = CurrentValue() = Matrix(5, bins);
 
 // Asigna los labels al objeto Matrix para cada columna
 m.SetLabel(0, "v");
 m.SetLabel(1, "vx dist");
 m.SetLabel(2, "vy dist");
 m.SetLabel(3, "vz dist");
 m.SetLabel(4, "|v| dist");
 for (int z=0;z<bins;++z)
 {
  m.Set(0, z, vrmin+z*(vrmax-vrmin)/double(bins-1));
  for (int q=0;q<3;++q) m.Set(q+1, z, 0.0);
  m.Set(4, z, 0.0);  
 }
 for (long int i=0;i<n;++i)
 {
  const Vector & v = atoms[i].Velocity();
  for (int q=0;q<3;++q)
  {
   int k = int(floor((bins-1)*((v[q]-vrmin)/(vrmax-vrmin))));
   m.Set(q+1, k, m.Get(q+1, k)+1);
  }
  int k = int(floor((bins-1)*((v.Module()-vrmin)/(vrmax-vrmin))));
  m.Set(4, k, m.Get(4, k)+1);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new VelDist(args); }
void destroy(Module * m) { delete m; }


