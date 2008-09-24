//
//
//

#ifndef __DLPOLYMODULE_H__
#define __DLPOLYMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class DlPolyFormat: public lpmd::CellFormat, public lpmd::Module
{
 public:
  //Metodos Generales
  DlPolyFormat(std::string args);
  virtual ~DlPolyFormat();
  void ShowHelp() const;
  std::string Keywords() const;

  //Metodos propios de modulo dlpoly
  void ReadHeader(std::istream & is) const;
  bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
  void WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *cells=NULL) const;
  void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
  long int GetInterval() const { return interval; }

 private:
  long int interval, *stepcnt;
  int level, pbkey;
  std::string ftype;
  bool rcell;
  double dt;
};

#endif

