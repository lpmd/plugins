//
//
//

#ifndef __CORDNUM_H__
#define __CORDNUM_H__

#include <map>

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class CordNum: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
    CordNum(double rcutoff, int nbins, double **cutoffs,int Ncutoffs);
    CordNum(std::string args);
    ~CordNum();

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
    std::map<std::string, double> rcut;
    int na;
    std::vector<std::string> satoms; //atomic symbols.
    bool do_average;
};

#endif

