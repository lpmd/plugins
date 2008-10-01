//
//
//

#include "cna.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/neighbor.h>
#include <lpmd/simulationcell.h>

#include <sstream>
#include <map>

using namespace lpmd;

class BondedPair
{
 public:
  BondedPair(unsigned int j0, unsigned int k0, unsigned int l0, double r0, const Vector & c): j(j0), k(k0), l(l0), r(r0), center(c) { }
  unsigned int j, k, l;
  double r;
  Vector center;
};

CommonNeighborAnalysis::CommonNeighborAnalysis(std::string args): Module("cna")
{
 m = NULL;
 AssignParameter("mode", "statistics");
 ProcessArguments(args);
 std::string m = GetString("mode");
 if (m == "full") mode = 0;
 else if (m == "statistics") mode = 1;
 rcut = GetDouble("rcut");
 start_step = GetInteger("start");
 end_step = GetInteger("end");
 interval = GetInteger("each");
 outputfile = GetString("output");
}

CommonNeighborAnalysis::~CommonNeighborAnalysis() { if (m != NULL) delete m; }

void CommonNeighborAnalysis::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = cna                                                      \n";
 std::cout << " Module Version     = 1.0                                                      \n";
 std::cout << " Support API lpmd   = 1.0.0                                                    \n";
 std::cout << " Problems Report to = gnm@gnm.cl                                               \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      Modulo utilizado para realizar Common Neighbor Analysis (CNA) a un       \n";
 std::cout << "      conjunto de configuraciones atomicas. Este metodo permite una determina- \n";
 std::cout << "      cion del tipo de estructura (FCC, BCC, HCP, etc) mas precisa que el uso  \n";
 std::cout << "      de la distribucion de pares g(r). Cada par de primeros vecinos es rotula-\n";
 std::cout << "      do segun 3 indices (j, k, l), donde:                                     \n";
 std::cout << "         j: numero de vecinos comunes a ambos atomos del par                   \n";
 std::cout << "         k: numero de enlaces que se pueden formar entre los j vecinos         \n";
 std::cout << "         l: longitud de la cadena continua mas larga formada por los k enlaces \n\n";
 std::cout << "      Como reconocer una estructura:                                           \n";
 std::cout << "      FCC perfecta: 100% pares de tipo 4-2-1                                   \n";
 std::cout << "      HCP perfecta: 50% pares tipo 4-2-1 y 50% pares tipo 4-2-2                \n";
 std::cout << "      BCC perfecta: 42.857% (3/7) de pares 4-4-4 y 57.143% (4/7) de pares 6-6-6\n\n";
 std::cout << "      El radio de corte (rcut) recomendado es el primer minimo de la g(r)      \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      rcut          : Especifica el radio maximo para el conteo de pares       \n";
 std::cout << "      mode          : Especifica el modo de salida (full/statistics)           \n";
 std::cout << "      output        : Fichero en el que se graba la salida                     \n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use cna                                                                       \n";
 std::cout << "     output cna.dat                                                            \n";
 std::cout << "     mode statistics                                                           \n";
 std::cout << "     rcut 2.0                                                                  \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al Modulo :                                                          \n";  
 std::cout << " property cna start=1 each=10 end=100                                          \n\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string CommonNeighborAnalysis::Keywords() const { return "rcut start end step output"; }

std::vector<long int> AddSegmentToChain(std::vector<long int> & cnm, std::map<long int, std::vector<long int> > & bonds, std::vector<long int> & chain, std::vector<long int> & longest)
{
 long int tail = chain[chain.size()-1];
 bool flag = false;
 for (unsigned int i=0;i<bonds[tail].size();++i)
 {
  bool inchain = false;
  for (unsigned int j=0;j<chain.size();++j) 
    if (chain[j] == bonds[tail][i]) 
    {
     inchain = true;
     break;
    }
  if ((chain.size() > 2) && (bonds[tail][i] == chain[0])) inchain = false; // exception to allow a closed chain
  if (! inchain) 
  {
   flag = true;
   std::vector <long int> newchain = chain, n;
   newchain.push_back(bonds[tail][i]);
   n = AddSegmentToChain(cnm, bonds, newchain, longest);
   if (n.size() > longest.size()) longest = n;
  }
 }
 if (! flag) return chain;
 std::vector<long int> emp;
 return emp;
}

void CommonNeighborAnalysis::Evaluate(SimulationCell & simcell, Potential & pot)
{
 const long n = simcell.Size();
 std::vector<BondedPair> data;
 std::list<Neighbor> * neighbormatrix = new std::list<Neighbor>[n];
 for (long i=0;i<n;++i) simcell.BuildNeighborList(i, neighbormatrix[i], true, rcut);
 for (long i=0;i<n;++i)
 {
  std::list<Neighbor> & nlist = neighbormatrix[i];
  for(std::list<Neighbor>::const_iterator it=nlist.begin();it!=nlist.end();++it)
  {
   const Neighbor & nn = *it;
   if (nn.r >= rcut) continue;  
   unsigned int cna_indices[3];
   long int j = nn.j->Index();
   std::list<Neighbor> & jlist = neighbormatrix[j];
   std::list<long int> cn;
   std::vector<long int> cnm;
   for (std::list<Neighbor>::const_iterator jt=jlist.begin();jt!=jlist.end();++jt)
   {
    if ((*jt).r >= rcut) continue;
    if (((*jt).j->Index() != i) && ((*jt).j->Index() != j)) cn.push_back((*jt).j->Index());
   }
   for (std::list<long int>::const_iterator qt=cn.begin();qt!=cn.end();++qt)
   {
    for (std::list<Neighbor>::const_iterator kt=nlist.begin();kt!=nlist.end();++kt)
     if (((*kt).j->Index() == (*qt)) || ((*kt).i->Index() == (*qt)))
     {
      cnm.push_back(*qt);
      break;
     }
   }
   cna_indices[0] = cnm.size(); // number of common neighbors
   // now count the number of bonds
   cna_indices[1] = 0;
   cna_indices[2] = 0;
   if (cna_indices[0] > 0) 
   {
    std::map<long int, std::vector<long int> > bonds;
    for (unsigned int p=0;p<cna_indices[0]-1;++p)
     for (unsigned int q=p+1;q<cna_indices[0];++q)
     {
      if (simcell.Distance(cnm[p], cnm[q]) < rcut)
      {
       bonds[cnm[p]].push_back(cnm[q]);
       bonds[cnm[q]].push_back(cnm[p]);
       cna_indices[1]++;
      }
     }
    // now try to connect the bonds in the longest chain ...
    std::vector<long int> longest;
    for (unsigned int p=0;p<cna_indices[0];++p)
    {
     std::vector<long int> tmp, n;
     tmp.push_back(cnm[p]);
     n = AddSegmentToChain(cnm, bonds, tmp, longest);
     if (n.size() > longest.size()) longest = n;
    }
    cna_indices[2] = longest.size()-1;
   }
   data.push_back(BondedPair(cna_indices[0], cna_indices[1], cna_indices[2], nn.r, simcell[i].Position()+nn.rij*0.5));
  }
 }
 delete [] neighbormatrix;
 unsigned long int npairs = data.size();
 if (m != NULL) delete m;
 if (mode == 0)
 {
  // Full mode
  m = new Matrix(7, npairs);
  // Asigna los labels al objeto Matrix para cada columna
  m->SetLabel(0, "j");
  m->SetLabel(1, "k");
  m->SetLabel(2, "l");
  m->SetLabel(3, "r");
  m->SetLabel(4, "x0");
  m->SetLabel(5, "y0");
  m->SetLabel(6, "z0");
  for (unsigned long int q=0;q<npairs;++q)
  {
   m->Set(0, q, data[q].j);
   m->Set(1, q, data[q].k);
   m->Set(2, q, data[q].l);
   m->Set(3, q, data[q].r);
   m->Set(4, q, data[q].center.Get(0));
   m->Set(5, q, data[q].center.Get(1));
   m->Set(6, q, data[q].center.Get(2));
  }
 } 
 else if (mode == 1)
 {
  // Statistics mode
  std::map<std::string, unsigned long int> stat;
  for (unsigned long int q=0;q<npairs;++q)
  {
   std::string s = ToString<int>(data[q].j)+"-"+ToString<int>(data[q].k)+"-"+ToString<int>(data[q].l);
   if (stat.count(s) == 0) stat[s] = 0;
   stat[s]++;
  }
  int nkeys = 0;
  for (std::map<std::string, unsigned long int>::const_iterator qt=stat.begin();qt!=stat.end();++qt) nkeys++;
  m = new Matrix(4, nkeys);
  m->SetLabel(0, "j");
  m->SetLabel(1, "k");
  m->SetLabel(2, "l");
  m->SetLabel(3, "percentage");
  int nk = 0;
  for (std::map<std::string, unsigned long int>::const_iterator qt=stat.begin();qt!=stat.end();++qt)
  {
   std::string s = (*qt).first;
   std::vector<std::string> splt = SplitTextLine(s, '-');
   m->Set(0, nk, atoi(splt[0].c_str()));
   m->Set(1, nk, atoi(splt[1].c_str()));
   m->Set(2, nk, atoi(splt[2].c_str()));
   m->Set(3, nk, 100.0*(double(stat[s])/double(npairs)));
   nk++;
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new CommonNeighborAnalysis(args); }
void destroy(Module * m) { delete m; }

