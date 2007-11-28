//
//
//

#ifndef __TESTMODULE_H__
#define __TESTMODULE_H__

#include <lpmd/plugin.h>

class TestModule: public lpmd::Module
{
 public:
  //
  TestModule(std::string args);
  virtual ~TestModule();

  void SetParameter(std::string key);
  void Show() const;
  std::string Keywords() const;
};

#endif



