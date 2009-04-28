//
//
//

#include "selectatoms.h"

#include <lpmd/md.h>
#include <lpmd/simulationcell.h>

#include <iostream>

using namespace lpmd;

SelectAtomsModifier::SelectAtomsModifier(std::string args): Module("selectatoms")
{
 AssignParameter("outside", "false");
 AssignParameter("from_index", "-1");
 AssignParameter("to_index", "-1");
 AssignParameter("xmin", "0.0");
 AssignParameter("xmax", "1.0");
 AssignParameter("ymin", "0.0");
 AssignParameter("ymax", "1.0");
 AssignParameter("zmin", "0.0");
 AssignParameter("zmax", "1.0");
 AssignParameter("mode", "extract");
 AssignParameter("vx", "0.0");
 AssignParameter("vy", "0.0");
 AssignParameter("vz", "0.0");
 // hasta aqui los valores por omision
 ProcessArguments(args);
 mode = GetString("mode");
 p0 = GetInteger("from_index");
 p1 = GetInteger("to_index");
 vmin[0] = GetDouble("xmin");
 vmax[0] = GetDouble("xmax");
 vmin[1] = GetDouble("ymin");
 vmax[1] = GetDouble("ymax");
 vmin[2] = GetDouble("zmin");
 vmax[2] = GetDouble("zmax");
 vel[0] = GetDouble("vx");
 vel[1] = GetDouble("vy");
 vel[2] = GetDouble("vz");
 outside = GetBool("outside");
}

SelectAtomsModifier::~SelectAtomsModifier() { }

void SelectAtomsModifier::SetParameter(std::string name)
{
 if (name == "index") 
 {
  AssignParameter("from_index", GetNextWord());
  AssignParameter("to_index", GetNextWord());
 }
 else if (name == "x")
 {
  AssignParameter("xmin", GetNextWord());
  AssignParameter("xmax", GetNextWord());
 }
 else if (name == "y")
 {
  AssignParameter("ymin", GetNextWord());
  AssignParameter("ymax", GetNextWord());
 }
 else if (name == "z")
 {
  AssignParameter("zmin", GetNextWord());
  AssignParameter("zmax", GetNextWord());
 }
 else if (name == "outside") AssignParameter("outside", "true");
 else if (name == "inside") AssignParameter("outside", "false");
 else Module::SetParameter(name);
}

void SelectAtomsModifier::Show(std::ostream & os) const
{
 Module::Show(os);
 if (p0 != -1)
 {
  os << "  from index = " << p0 << '\n';
  os << "  to index   = " << p1 << '\n';
 }
 else
 {
  if (mode == "extract")
  {
   if (outside) os << "Extracting atoms NOT inside range: " << '\n';
   else os << "Extracting atoms inside range: " << '\n';
  }
  if (mode == "setvelocity")
  {
   if (outside) os << "Setting velocity to atoms NOT inside range: " << '\n';
   else os << "Setting velocity to atoms inside range: " << '\n';
  }
  os << "  x between " << vmin[0] << " and " << vmax[0] << '\n';
  os << "  y between " << vmin[1] << " and " << vmax[1] << '\n';
  os << "  z between " << vmin[2] << " and " << vmax[2] << '\n';
  if (mode == "setvelocity") os << "Velocity set to " << vel << '\n';
 }
}

void SelectAtomsModifier::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << '\n';
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use selectatoms                                                               \n";
 std::cout << "     index 50 75                                                               \n";
 std::cout << " enduse                                                                        \n\n";
 std::cout << " use selectatoms                                                               \n";
 std::cout << "     x 10.0 15.0                                                               \n";
 std::cout << "     y 12.0 13.0                                                               \n";
 std::cout << "     z 2.0 8.0                                                                 \n";
 std::cout << " enduse                                                                        \n\n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " apply selectatoms                                                             \n\n";
}

void SelectAtomsModifier::Apply(SimulationCell & sc)
{
 const long n = sc.size();
 if (mode == "extract")
 {
  ParticleSet tmp;
  tmp.clear();
  if (p0 != -1)
  {
   for (long i=p0;i<=p1;++i) tmp.Create(new Atom(sc[i]));
  }
  else
  {
   for (long i=0;i<n;++i)
   {
    const Vector & pos = sc[i].Position();
    bool select_this = false;
    if ((pos[0] >= vmin[0]) && (pos[0] <= vmax[0]))
    {
     if ((pos[1] >= vmin[1]) && (pos[1] <= vmax[1]))
     {
      if ((pos[2] >= vmin[2]) && (pos[2] <= vmax[2])) select_this = true;
     }
    }
    if (outside) select_this = (!select_this);
    if (select_this) tmp.Create(new Atom(sc[i]));
   }
  }
  sc.clear();
  for (unsigned long i=0;i<tmp.size();++i) sc.Create(new Atom(tmp[i]));
  sc.NumEspec();
  sc.AssignIndex();
 }
 else if (mode == "setvelocity")
 {
  if (p0 != -1) 
  {
   for (long i=p0;i<=p1;++i) sc.SetVelocity(i, vel);
  } 
  else
  {
   for (long i=0;i<n;++i)
   {
    const Vector & pos = sc[i].Position();
    bool select_this = false;
    if ((pos[0] >= vmin[0]) && (pos[0] <= vmax[0]))
    {
     if ((pos[1] >= vmin[1]) && (pos[1] <= vmax[1]))
     {
      if ((pos[2] >= vmin[2]) && (pos[2] <= vmax[2])) select_this = true;
     }
    }
    if (outside) select_this = (!select_this);
    if (select_this) sc.SetVelocity(i, vel);
   }
  }
 }
}

void SelectAtomsModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new SelectAtomsModifier(args); }
void destroy(Module * m) { delete m; }

