//
//
//

#ifndef __PDBMODULE_H__
#define __PDBMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class PDBFormat: public lpmd::CellFormat, public lpmd::Module
{
 public:
   //Metodos Generales
   PDBFormat(std::string args);
   virtual ~PDBFormat();
   void ShowHelp() const;

   //Metodos Propios de modulo mol2
   void WriteHeader(std::ostream & os, std::vector<lpmd::SimulationCell> *) const;
   void ReadHeader(std::istream & is) const;
   void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
   bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
   
   long int GetInterval() const {return interval;}

 private:
   long int interval;
};

#endif

