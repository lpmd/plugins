//
//
//

#ifndef __ANGDIST_H__
#define __ANGDIST_H__

#include <lpmd/scalartable.h>
#include <lpmd/instantproperty.h>
#include <lpmd/plugin.h>

#include <map>

class AngDist: public lpmd::ScalarTable, public lpmd::InstantProperty, public lpmd::Module
{
 public:
    //Metodos Generales
    AngDist(std::string args);
    ~AngDist();

    void SetParameter(std::string name);
    void Show(std::ostream & os) const;
    void ShowHelp() const;
    std::string Keywords() const;

    //Metodos Propios del modulo angdist
    const lpmd::Matrix & Value() const { return *m; }
    void Evaluate(lpmd::SimulationCell & simcell, lpmd::Potential & pot);

 private:
    lpmd::Matrix * m;
    int nb;
    std::map<std::string, double> rcut;
    int na;
    bool do_average;
    std::vector<std::string> satoms; //simbolos de los atomos.
};

#endif
