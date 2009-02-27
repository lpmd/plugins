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

  //Metodos propios de modulo xyz
  void ReadHeader(std::istream & is) const;
  bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
  void WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *) const;
  void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
  long int GetInterval() const { return interval; }

 private:
  long int * linecounter;
  long int interval;
  long int level;
  bool rcell;
  std::string coords;
  std::string inside;
  std::string external;
};

#endif

