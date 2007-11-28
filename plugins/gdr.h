//
//
//

#ifndef __GDR_H__
#define __GDR_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class Gdr: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
    Gdr(std::string args);
    Gdr(double rcutoff, int nbins);
    ~Gdr();

    const lpmd::Matrix & Value() const { return *m; }

    // From Module
    void SetParameter(std::string name);
    void Show() const;
    void ShowHelp() const;
    std::string Keywords() const;

    // From InstantProperty
    void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    double rcut;
    int nb;
    bool do_average;
};

#endif

