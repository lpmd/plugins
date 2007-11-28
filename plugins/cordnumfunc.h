//
//
//

#ifndef __CORDNUMFUNC_H__
#define __CORDNUMFUNC_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class CordNumFunc: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
    CordNumFunc(double rcutoff, int nbins, double cutoff,int Ncutoffs);
    CordNumFunc(std::string args);
    ~CordNumFunc();

    // From Module
    void SetParameter(std::string name);
    void Show() const;
    void ShowHelp() const;
    std::string Keywords() const;

    const lpmd::Matrix & Value() const { return *m; }
    void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int nb;
    double cut;
    int na;
    std::vector<std::string> satoms; //atomic symbols.
    bool do_average;
};

#endif

