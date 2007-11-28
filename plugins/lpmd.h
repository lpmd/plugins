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
   //
   LPMDFormat(std::string args);
   virtual ~LPMDFormat();

   void SetParameter(std::string name);
   void Show() const;
   void ShowHelp() const;
   std::string Keywords() const;

   // Overriden from CellReader
   void ReadHeader(std::istream & is) const;
   void ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;

   // Overriden from CellWriter
   void WriteHeader(std::ostream & os) const;
   void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
   
   long int GetInterval() const {return interval;}

 private:
   long int interval;
   long int level;
};

#endif

