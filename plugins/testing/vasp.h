//
//
//

#ifndef __VASPMODULE_H__
#define __VASPMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class VaspFormat: public lpmd::CellFormat, public lpmd::Plugin
{
 public:
  //Metodos Generales
  VaspFormat(std::string args);
  virtual ~VaspFormat();
  void ShowHelp() const;

  //Metodos propios de modulo dlpoly
  void ReadHeader(std::istream & is) const;
  bool ReadCell(std::istream & is, lpmd::Configuration & con) const;
  void WriteHeader(std::ostream & os, lpmd::SimulationHistory * sh) const;
  void WriteCell(std::ostream & os, lpmd::Configuration & con) const;
  long int GetInterval() const { return interval; }

 public:
  long int interval;
  int level;
  bool rcell;
  std::string speclist;
  lpmd::Array<std::string> satoms; //simbolos atomicos
  std::string tp; // tipo de celda Direct/Cartesian
};

#endif

