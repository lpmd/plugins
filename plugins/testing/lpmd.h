//
//
//

#ifndef __LPMDMODULE_H__
#define __LPMDMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>
#include <lpmd/simulation.h>
#include <lpmd/array.h>
#include <string>

class LPMDFormat: public lpmd::CellFormat, public lpmd::Plugin
{
 public:
   //Metodos Generales
   LPMDFormat(std::string args);
   virtual ~LPMDFormat();
   void ShowHelp() const;

   //Metodos Propios de modulo lpmd
   void ReadHeader(std::istream & is) const;
   bool ReadCell(std::istream & is, lpmd::Configuration & con) const;
   void WriteHeader(std::ostream & os, lpmd::SimulationHistory *) const;
   void WriteCell(std::ostream & os, lpmd::Configuration & con) const;
   
   long int GetInterval() const {return interval;}

 private:
   long int * linecounter;
   long int interval;
   mutable long int level;
   lpmd::Array<std::string> extra;
   bool rcell;
   mutable lpmd::Array<std::string> hdr;
};

#endif

