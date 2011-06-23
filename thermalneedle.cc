//
//
//

#include "thermalneedle.h"

#include <lpmd/simulationcell.h>
#include <lpmd/util.h>
#include <lpmd/md.h>

#include <iostream>

using namespace lpmd;

ThermalNeedleModifier::ThermalNeedleModifier(std::string args): Module("tempscaling")
{
 AssignParameter("temp", "300.0");
 AssignParameter("center", "<0.5,0.5,0.5> ");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 temperature = GetDouble("temp");
 center = GetVector("center");
 radius = GetDouble("radius");
 start = GetInteger("start");
 end = GetInteger("end");
 each = GetInteger("each");
}

ThermalNeedleModifier::~ThermalNeedleModifier() { }

void ThermalNeedleModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo es utilizado para escalar la temperatura de un grupo de        \n";
 std::cout << " particulas, en una zona esferica ubicada en 'center' y de radio 'radius'.     \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      temp          : Temperatura a aplicar en la zona.                        \n";
 std::cout << "      center        : Centro donde se ubica la esfera, debe ser entregado en   \n";
 std::cout << "                      formato de vector '<X,Y,Z>'.                             \n";
 std::cout << "      radius        : Radio de la esfera en donde se aplica la temperatura.    \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use thermalneedle                                                             \n";
 std::cout << "     temp 300                                                                  \n";
 std::cout << "     center <14,15.5,23.6>                                                     \n";
 std::cout << "     radius 5                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo                                                            \n";
 std::cout << " apply thermalneedle start=1000 each=10 end=1000                             \n\n";
 std::cout << "      De esta forma aplicamos el termostato entre 100 y 1000 cada 10 steps.    \n";
}

void ThermalNeedleModifier::Apply(SimulationCell & sc)
{
 //sc.SetTemperature(fromtemp);
}

void ThermalNeedleModifier::Apply(MD & md)
{
 SimulationCell & sc = md.GetCell();
 std::vector<int> idatoms;
 //Nuevo Set de particulas dentro de la region
 SimulationCell parts;              // FIXME 0.5
 // long int i=0;
 for(unsigned long i=0;i<=sc.size();++i)
 {
  const Vector & pos = sc[i].Position();
  double r = sc.Displacement(center, pos).Module();
  if ( r < radius ) 
  {
   parts.Create(new Atom(sc[i]));   // FIXME 0.5 
   idatoms.push_back(i);
  }
 }
 //Reasignamos nueva temperatura a las particulas dentro de esa region.
 parts.SetTemperature(temperature); 
 //Reemplazamos las particulas de la lista idatoms
 //con las de parts.
 for(unsigned long int i=0 ;i<parts.size() ; ++i)
 {
  int id = idatoms[i];
  sc[id] = parts[i];    // FIXME 0.5
 }
 std::cerr << "-> Apply ThemralNeedle to T = " << temperature << '\n';  
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new ThermalNeedleModifier(args); }
void destroy(Module * m) { delete m; }
