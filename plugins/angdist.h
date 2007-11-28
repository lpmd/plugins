//
//
//

#ifndef __ANGDIST_H__
#define __ANGDIST_H__

#include <map>

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

class AngDist: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
    AngDist(double rcutoff, int nbins, int Ncutoffs);
    AngDist(std::string args);
    ~AngDist();

    void SetParameter(std::string name);
    void Show() const;
    std::string Keywords() const;

    const lpmd::Matrix & Value() const { return *m; }
    void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int nb;
    std::map<std::string, double> rcut;
    //double **cut;
    int na;
    bool do_average;
    std::vector<std::string> satoms; //atomic symbols.
};

#endif

