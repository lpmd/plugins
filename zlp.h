//
//
//

#ifndef __ZLPMODULE_H__
#define __ZLPMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/simulationhistory.h>
#include <lpmd/plugin.h>

using namespace lpmd;

class ZLPFormat: public CellFormat, public Plugin
{
 public:
   //Metodos Generales
   ZLPFormat(std::string args);
   virtual ~ZLPFormat();
   void ShowHelp() const;

   //Metodos Propios de modulo zlp
   void ReadHeader(std::istream & is) const;
   bool ReadCell(std::istream & is, Configuration & conf) const;
   void WriteHeader(std::ostream & os, SimulationHistory * sh) const;
   void WriteCell(std::ostream & os, Configuration & conf) const;
   
   long int GetInterval() const { return interval; }

 private:
   void * zstr;            // z_stream structure, used with zlib 
   unsigned char * inbuf, * outbuf;
   long int interval;
   int level, blocksize, complev, *lastop;
   bool rcell;
};

#endif

