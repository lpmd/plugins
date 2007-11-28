//
//
//

#include <iostream>

#include "test.h"

using namespace lpmd;

TestModule::TestModule(std::string args): Module("test") 
{ 
 ProcessArguments(args); 
}

TestModule::~TestModule() { }

void TestModule::SetParameter(std::string key)
{
 if (key == "arg1") AssignParameter("arg1", GetNextWord());
 if (key == "arg2") AssignParameter("arg2", GetNextWord());
 if (key == "magicnumber") AssignParameter("magicnumber", GetNextWord());
}

void TestModule::Show() const
{
 Module::Show();
 std::cout << "   arg1 = " << GetString("arg1") << '\n';
 std::cout << "   arg2 = " << GetString("arg2") << '\n';
 std::cout << "   magicnumber = " << GetString("magicnumber") << '\n';
}

std::string TestModule::Keywords() const
{
 return "arg1 arg2 magicnumber";
}

// This is included so that the module can be loaded dynamically
Module * create(std::string args) { return new TestModule(args); }
void destroy(Module * m) { delete m; }


