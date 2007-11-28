//
//
//

#ifndef __POVRAY_H__
#define __POVRAY_H__

#include <lpmd/visualizer.h>
#include <lpmd/plugin.h>

class POVRAY: public lpmd::Visualizer, public lpmd::Module
{
 public:
    POVRAY(std::string args);
    ~POVRAY();

    // From Module
    void SetParameter(std::string name);
    void Show() const;
    void ShowHelp() const;
    std::string Keywords() const;

    // From Visualizer
    void Apply(const lpmd::MD & md);

 private:
    std::string header;
    std::string direct;
    int counter;
    int level;
};

#endif

