//
//
//

#ifndef __LPMDMODULE_H__
#define __LPMDMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class LPMDFormat: public lpmd::CellFormat, public lpmd::Module
{
 public:
   //Metodos Generales
   LPMDFormat(std::string args);
   virtual ~LPMDFormat();
   void ShowHelp() const;
   std::string Keywords() const;

   //Metodos Propios de modulo lpmd
   void ReadHeader(std::istream & is) const;
   bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
   void WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *cell=NULL) const;
   void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
   
   long int GetInterval() const {return interval;}

 private:
   long int * linecounter;
   long int interval;
   long int level;
   bool rcell;
};

#endif

