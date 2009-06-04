//
//
//

#ifndef __SPHERE_FILTER_H__
#define __SPHERE_FILTER_H__

#include <lpmd/systemfilter.h>
#include <lpmd/atomselector.h>
#include <lpmd/plugin.h>

class SphereFilter: public lpmd::SystemFilter, public lpmd::Module
{
 public:
  SphereFilter(std::string args);
  ~SphereFilter();

  void ShowHelp() const;

  lpmd::Selector<lpmd::BasicParticleSet> & CreateSelector();

  private:
   double radius;
   lpmd::Vector center;
   lpmd::Selector<lpmd::BasicParticleSet> * selector;
};

#endif

