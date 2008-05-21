//
//
//

#ifndef __ZLPMODULE_H__
#define __ZLPMODULE_H__

#include <lpmd/cellformat.h>
#include <lpmd/plugin.h>

class ZLPFormat: public lpmd::CellFormat, public lpmd::Module
{
 public:
   //Metodos Generales
   ZLPFormat(std::string args);
   virtual ~ZLPFormat();
   void ShowHelp() const;
   std::string Keywords() const;

   //Metodos Propios de modulo zlp
   void ReadHeader(std::istream & is) const;
   bool ReadCell(std::istream & is, lpmd::SimulationCell & sc) const;
   void WriteHeader(std::ostream & os) const;
   void WriteCell(std::ostream & os, lpmd::SimulationCell & sc) const;
   
   long int GetInterval() const {return interval;}

 private:
   void * zstr;            // z_stream structure, used with zlib 
   unsigned char * inbuf, * outbuf;
   long int interval;
   int level, blocksize, complev, *lastop;
};

#endif

