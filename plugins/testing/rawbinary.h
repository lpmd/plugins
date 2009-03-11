//
//
//

#ifndef __RAWBINMODULE_H__
#define __RAWBINMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>
#include <string.h>

class RawBinFormat: public lpmd::CellFormat, public lpmd::Module
{
 public:
  //Metodos Generales
  RawBinFormat(std::string args);
  virtual ~RawBinFormat();
  void ShowHelp() const;
  std::string Keywords() const;

  void ReadHeader(std::istream & is) const;
  bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
  void WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *cell=NULL) const;
  void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
  long int GetInterval() const { return interval; }

 private:
  long int interval;
  long int level;
  bool rcell;
};

#endif

