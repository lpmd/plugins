//
//
//

#include "cna.h"

#include <lpmd/matrix.h>
#include <lpmd/util.h>
#include <lpmd/atompair.h>
#include <lpmd/configuration.h>

#include <sstream>
#include <map>

using namespace lpmd;

class BondedPair: public IndexTrio
{
 public:
  BondedPair(unsigned int j0, unsigned int k0, unsigned int l0, unsigned long int ii, unsigned long int jj, double r0, const Vector & c): IndexTrio(j0, k0, l0), ati(ii), atj(jj), r(r0), center(c) { }
  unsigned int ati, atj; // index of atoms involved in the pair
  double r;
  Vector center;
};

CommonNeighborAnalysis::CommonNeighborAnalysis(std::string args): Plugin("cna", "2.0")
{
 ParamList & params = (*this);
 //
 DefineKeyword("rcut");
 DefineKeyword("reference");
 DefineKeyword("mode", "statistics");
 DefineKeyword("species", "all");
 DefineKeyword("filterby", "none");
 DefineKeyword("start", "0");
 DefineKeyword("end", "-1");
 DefineKeyword("each", "1");
 DefineKeyword("output");
 DefineKeyword("filterout", "all");
 //
 ProcessArguments(args);
 std::string m = params["mode"];
 if (m == "full") mode = 0;
 else if (m == "statistics") mode = 1;
 else if (m == "defects") mode = 2;
 rcut = double(params["rcut"]);
 std::string rfs = params["reference"];
 if (rfs == "fcc") refmap[IndexTrio(4, 2, 1)] = 1;
 else if (rfs == "bcc")
 {
  refmap[IndexTrio(4, 4, 4)] = 1;
  refmap[IndexTrio(6, 6, 6)] = 1;
 }
 else if (rfs == "hcp")
 {
  refmap[IndexTrio(4, 2, 1)] = 1;
  refmap[IndexTrio(4, 2, 2)] = 1;
 }
 std::string species = params["species"];
 spc1 = spc2 = -1;
 if (species != "all")
 {
  Array<std::string> tmp = StringSplit(species, '-');
  spc1 = ElemNum(tmp[0]); 
  spc2 = ElemNum(tmp[1]);
 }
 start = int(params["start"]);
 end = int(params["end"]);
 each = int(params["each"]);
 OutputFile() = params["output"];
 filterout = params["filterout"];
}

CommonNeighborAnalysis::~CommonNeighborAnalysis() { }

void CommonNeighborAnalysis::ShowHelp() const
{
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
 std::cout << "      mode          : Especifica el modo de salida (full/statistics/defects)   \n";
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
}

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

void CommonNeighborAnalysis::Evaluate(Configuration & conf, Potential & pot)
{
 assert(&pot != 0);//icc869
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 const long n = atoms.Size();
 std::vector<BondedPair> data;
 NeighborList * neighbormatrix = new NeighborList[n];
 DebugStream() << "-> Computing Common Neighbor Analysis (CNA) over " << n << " atoms\n";
 DebugStream() << "-> Building neighbor lists\n";
 for (long i=0;i<n;++i)
 {
  neighbormatrix[i] = conf.Neighbors(i, true, rcut);
 }
 // Construct an "index table" so we don't have to depend on Atom::Index()
 std::map<BasicAtom *, long int> indices;
 for (long int i=0;i<atoms.Size();++i) indices[&atoms[i]] = i;
 // 

 DebugStream() << "-> Searching for pairs, species = " << (*this)["species"] << '\n'; 
 for (long i=0;i<n;++i)
 {
  NeighborList & nlist = neighbormatrix[i];
  for (long int it=0;it<nlist.Size();++it)
  {
   const AtomPair & nn = nlist[it];
   if (spc1 != -1)
   {
    bool rightpair = (atoms[i].Z() == spc1) && (nn.j->Z() == spc2);
    rightpair = (rightpair || ((atoms[i].Z() == spc2) && (nn.j->Z() == spc1)));
    if (!rightpair) continue;
   }
   if (nn.r2 >= rcut*rcut) continue;  
   unsigned int cna_indices[3];
   long int j = indices[nn.j];
   NeighborList & jlist = neighbormatrix[j];
   std::list<long int> cn;
   std::vector<long int> cnm;
   for (long int jt=0;jt<jlist.Size();++jt)
   {
    const AtomPair & jnn = jlist[jt];
    if (jnn.r2 >= rcut*rcut) continue;
    if ((indices[jnn.j] != i) && (indices[jnn.j] != j)) cn.push_back(indices[jnn.j]);
   }
   for (std::list<long int>::const_iterator qt=cn.begin();qt!=cn.end();++qt)
   {
    for (long int kt=0;kt<nlist.Size();++kt)
     if ((indices[nlist[kt].j] == (*qt)) || (indices[nlist[kt].i] == (*qt)))
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
      if (cell.Displacement(atoms[cnm[p]].Position(), atoms[cnm[q]].Position()).Module() < rcut)
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
     std::vector<long int> tmp, nt;
     tmp.push_back(cnm[p]);
     nt = AddSegmentToChain(cnm, bonds, tmp, longest);
     if (nt.size() > longest.size()) longest = nt;
    }
    cna_indices[2] = longest.size()-1;
   }
   data.push_back(BondedPair(cna_indices[0], cna_indices[1], cna_indices[2], i, indices[nn.j], sqrt(nn.r2), atoms[i].Position()+nn.rij*0.5));
  }
 }
 delete [] neighbormatrix;
 unsigned long int npairs = data.size();

 Matrix * m = &(CurrentValue());
 if (mode == 0)
 {
  // Full mode
  *m = Matrix(7, npairs);
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
   m->Set(4, q, data[q].center[0]);
   m->Set(5, q, data[q].center[1]);
   m->Set(6, q, data[q].center[2]);
  }
 } 
 else if (mode == 1)
 {
  // Statistics mode
  std::map<std::string, unsigned long int> stat;
  std::map<std::string, double> stat_rav;
  std::map<std::string, double> stat_r2av;
  for (unsigned long int q=0;q<npairs;++q)
  {
   std::string s = ToString<int>(data[q].j)+"-"+ToString<int>(data[q].k)+"-"+ToString<int>(data[q].l);
   if (stat.count(s) == 0) 
   {
    stat[s] = 0;
    stat_rav[s] = 0.0;
    stat_r2av[s] = 0.0;
   }
   stat[s]++;
   stat_rav[s] += data[q].r;
   stat_r2av[s] += pow(data[q].r, 2.0);
  }
  int nkeys = 0;
  for (std::map<std::string, unsigned long int>::const_iterator qt=stat.begin();qt!=stat.end();++qt) nkeys++;
  *m = Matrix(6, nkeys);
  m->SetLabel(0, "j");
  m->SetLabel(1, "k");
  m->SetLabel(2, "l");
  m->SetLabel(3, "percentage");
  m->SetLabel(4, "R_average");
  m->SetLabel(5, "R_stddev");
  int nk = 0;
  for (std::map<std::string, unsigned long int>::const_iterator qt=stat.begin();qt!=stat.end();++qt)
  {
   std::string s = (*qt).first;
   double ns = double(stat[s]);
   Array<std::string> splt = StringSplit(s, '-');
   m->Set(0, nk, atoi(splt[0].c_str()));
   m->Set(1, nk, atoi(splt[1].c_str()));
   m->Set(2, nk, atoi(splt[2].c_str()));
   m->Set(3, nk, 100.0*(ns/double(npairs)));
   m->Set(4, nk, stat_rav[s]/ns);
   m->Set(5, nk, sqrt(stat_r2av[s]/ns-pow(stat_rav[s]/ns, 2.0)));
   nk++;
  }
  if ((stat.size() == 1) && ((*(stat.begin())).first == "4-2-1")) 
     DebugStream() << "-> All pairs 4-2-1, this is definitely FCC!\n";
  else if (stat.size() == 2)
  {
   Array<std::string> lista;
   std::map<std::string, unsigned long int>::iterator it = stat.begin();
   lista.Append((*it).first);
   lista.Append((*(++it)).first);
   if ((lista.Find("4-2-1") >= 0) && (lista.Find("4-2-2") >= 0)) DebugStream() << "-> Found 4-2-1 and 4-2-2 pairs, looks like HCP\n";
   if ((lista.Find("6-6-6") >= 0) && (lista.Find("4-4-4") >= 0)) DebugStream() << "-> Found 6-6-6 and 4-4-4 pairs, looks like BCC\n";
  }
  else DebugStream() << "-> Found " << stat.size() << " kinds of pairs, cannot tell the structure\n";
 }
 else if (mode == 2)
 {
  // Defects mode
  unsigned long int * regcnt = new unsigned long int[atoms.Size()];
  unsigned long int * defcnt = new unsigned long int[atoms.Size()];
  int nfil = 0;
  lpmd::Atom at;
  for (long int q=0;q<atoms.Size();++q) 
  {
   at = atoms[q];
   regcnt[q] = defcnt[q] = 0;
   if ((atoms.Have(at, Tag(filterout)) && (atoms.GetTag(at, Tag(filterout))
       == "true")) || (filterout=="all")) nfil++;
  }
  for (unsigned long int q=0;q<npairs;++q)
  {
   if (refmap.count(IndexTrio(data[q].j, data[q].k, data[q].l)) > 0) 
   {
    regcnt[data[q].ati]++;
    regcnt[data[q].atj]++;
   }
   else
   {
    defcnt[data[q].ati]++;
    defcnt[data[q].atj]++;
   }
  }
  *m = Matrix(6, nfil);
  m->SetLabel(0, "x");
  m->SetLabel(1, "y");
  m->SetLabel(2, "z");
  m->SetLabel(3, "reference");
  m->SetLabel(4, "defects");
  m->SetLabel(5, "%defect");
  for (long int q=0;q<atoms.Size();++q)
  {
   at = atoms[q];
   if ((atoms.Have(at, Tag(filterout)) && (atoms.GetTag(at, Tag(filterout))
      == "true")) || (filterout == "all"))
   {
    const Vector & pos = atoms[q].Position();
    for (int pp=0;pp<3;++pp) m->Set(pp, q, pos[pp]);
    m->Set(3, q, (double)regcnt[q]);
    m->Set(4, q, (double)defcnt[q]);
    m->Set(5, q, 100.0*defcnt[q]/(regcnt[q]+defcnt[q]));
   }
  }
  delete [] regcnt;
  delete [] defcnt;
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new CommonNeighborAnalysis(args); }
void destroy(Plugin * m) { delete m; }

