//
//
//

#ifndef __XYZMODULE_H__
#define __XYZMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class XYZFormat: public lpmd::CellFormat, public lpmd::Module
{
 public:
  //Metodos Generales
  XYZFormat(std::string args);
  virtual ~XYZFormat();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos propios de modulo xyz
  void ReadHeader(std::istream & is) const;
  bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
  void WriteHeader(std::ostream & os) const;
  void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
  long int GetInterval() const { return interval; }

 private:
  long int * linecounter;
  long int interval;
  long int level;
  std::string coords;
  std::string inside;
};

#endif

