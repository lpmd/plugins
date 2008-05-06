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
 // hasta aqui los valores por omision
 ProcessArguments(args);
 p0 = GetInteger("from_index");
 p1 = GetInteger("to_index");
 vmin.Set(0, GetDouble("xmin"));
 vmax.Set(0, GetDouble("xmax"));
 vmin.Set(1, GetDouble("ymin"));
 vmax.Set(1, GetDouble("ymax"));
 vmin.Set(2, GetDouble("zmin"));
 vmax.Set(2, GetDouble("zmax"));
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
  if (outside) os << "Selecting atoms NOT inside range: " << '\n';
  else os << "Selecting atoms inside range: " << '\n';
  os << "  x between " << vmin.Get(0) << " and " << vmax.Get(0) << '\n';
  os << "  y between " << vmin.Get(1) << " and " << vmax.Get(1) << '\n';
  os << "  z between " << vmin.Get(2) << " and " << vmax.Get(2) << '\n';
 }
}

void SelectAtomsModifier::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = selectatoms                                              \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string SelectAtomsModifier::Keywords() const
{
 return "index x y z";
}

void SelectAtomsModifier::Apply(SimulationCell & sc)
{
 Particles tmp;
 tmp.Clear();
 const long n = sc.Size();
 if (p0 != -1)
 {
  for (long i=p0;i<=p1;++i) tmp.AppendAtom(sc[i]);
 }
 else
 {
  for (long i=0;i<n;++i)
  {
   const Vector & pos = sc[i].Position();
   bool select_this = false;
   if ((pos.Get(0) >= vmin.Get(0)) && (pos.Get(0) <= vmax.Get(0)))
   {
    if ((pos.Get(1) >= vmin.Get(1)) && (pos.Get(1) <= vmax.Get(1)))
    {
     if ((pos.Get(2) >= vmin.Get(2)) && (pos.Get(2) <= vmax.Get(2))) select_this = true;
    }
   }
   if (outside) select_this = (!select_this);
   if (select_this) tmp.AppendAtom(sc[i]);
  }
 }
 sc.Clear();
 for (long i=0;i<tmp.Size();++i) sc.AppendAtom(tmp[i]);
 sc.NumEspec();
 sc.AssignIndex();
}

void SelectAtomsModifier::Apply(MD & md) { Apply(md.GetCell()); }

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new SelectAtomsModifier(args); }
void destroy(Module * m) { delete m; }

