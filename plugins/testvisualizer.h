//
//
//

#ifndef __TESTVISUALIZER_H__
#define __TESTVISUALIZER_H__

#include <lpmd/visualizer.h>
#include <lpmd/plugin.h>

class TestVisualizer: public lpmd::Visualizer, public lpmd::Module
{
 public:

   TestVisualizer(std::string args);
   ~TestVisualizer();

   //
   void SetParameter(std::string name);
   void Show() const;
   std::string Keywords() const;

   void Apply(const lpmd::MD & md);

};

#endif



