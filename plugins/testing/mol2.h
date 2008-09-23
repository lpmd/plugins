//
//
//

#ifndef __MOL2MODULE_H__
#define __MOL2MODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class Mol2Format: public lpmd::CellFormat, public lpmd::Module
{
 public:
   //Metodos Generales
   Mol2Format(std::string args);
   virtual ~Mol2Format();
   void ShowHelp() const;
   std::string Keywords() const;

   //Metodos Propios de modulo mol2
   void WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *cell=NULL) const;
   void ReadHeader(std::istream & is) const;
   void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
   bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
   
   long int GetInterval() const {return interval;}

 private:
   long int interval;
   bool rcell;
};

#endif

