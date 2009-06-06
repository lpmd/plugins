//
//
//

#include "lc2.h"

using namespace lpmd;

static double box[3] = { 1.0, 1.0, 1.0 };
static Configuration * isc = NULL;

double LinkedCellCellManager2::DistanceFunction(double dx, double dy, double dz)
{
 double dr[3];
 dr[0] = dx;
 dr[1] = dy;
 dr[2] = dz;
 for (int q=0;q<3;++q) 
 {
  double ll = box[q];
  if (dr[q] >= 0.5*ll) dr[q] -= ll;
  else if (dr[q] < -0.5*ll) dr[q] += ll;
 }
 return (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
}

LinkedCellCellManager2::LinkedCellCellManager2(std::string args): Plugin("lc2", "1.0")
{ 
 ParamList & params = (*this);
 AssignParameter("cutoff", "7.0");
 // 
 ProcessArguments(args);
 cutoff = params["cutoff"];
 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);

 LC_Init(&lc, nx, ny, nz); 
 pos = NULL;
}

LinkedCellCellManager2::~LinkedCellCellManager2() 
{ 
 LC_Destroy(&lc); 
 if (pos != NULL) free(pos);
}

void LinkedCellCellManager2::Show(std::ostream & os) const
{
 Module::Show(os);
}

void LinkedCellCellManager2::Reset() { }

void LinkedCellCellManager2::UpdateCell(Configuration & sc) 
{ 
 //std::cerr << "DEBUG in lc2 UpdateCell\n";
 BasicCell & cell = sc.Cell();
 BasicParticleSet & atoms = sc.Atoms();
 isc = &sc;
 for (int q=0;q<3;++q) box[q] = cell[q].Module();
 double (*dist2)(double,double,double) = (&LinkedCellCellManager2::DistanceFunction);
 LC_SetDistanceSqrFunction(&lc, dist2);
 if (pos == NULL) pos = (double *)(malloc(sizeof(double)*(3*atoms.Size())));
 for (long int i=0;i<atoms.Size();++i)
 {
  const Vector & rpos = atoms[i].Position();
  //std::cerr << "DEBUG rpos = " << rpos << '\n';
  Vector vpos = cell.Fractional(rpos);
  //std::cerr << "DEBUG vpos = " << vpos << '\n';
  /*
  assert(vpos[0] >= 0.0);
  assert(vpos[0] <= 1.0);
  assert(vpos[1] >= 0.0);
  assert(vpos[1] <= 1.0);
  assert(vpos[2] >= 0.0);
  assert(vpos[2] <= 1.0);
  */
  vpos = cell.ScaleByCell(vpos);
  for (int q=0;q<3;++q) pos[i*3+q] = vpos[q];
 }
 //std::cerr << "DEBUG in lc2 before LC_FillCells\n";
 LC_FillCells(&lc, atoms.Size(), pos, box);
 //std::cerr << "DEBUG in lc2 after LC_FillCells\n";
}

void LinkedCellCellManager2::BuildNeighborList(Configuration & sc, long i, NeighborList & nlist, bool full, double rcut) 
{ 
 //std::cerr << "DEBUG in lc2 BuildNeighborList\n";
 nlist.Clear();
 BasicParticleSet & atoms = sc.Atoms();
 BasicCell & cell = sc.Cell();
 NeighborInfo * p = 0, * neighbors = LC_GetNeighbors(&lc, i, rcut, pos, box, (full ? 1 : 0));
 AtomPair nn;
 nn.i = &atoms[i];
 for (p=neighbors;(p->j)>=0;++p)
 { 
  if (p->j == -1) break;
  nn.j = &atoms[p->j];
  //nn.rij = Vector(p->dr[0], p->dr[1], p->dr[2]);
  nn.rij = cell.Displacement(nn.i->Position(), atoms[p->j].Position());
  //std::cerr << "DEBUG " << Vector(p->dr[0], p->dr[1], p->dr[2]) << "   " << nn.rij << '\n';
  //nn.r = sqrt(p->r2);
  nn.r = nn.rij.Module();
  nlist.Append(nn);
 }
 free(neighbors);
 //std::cerr << "DEBUG in lc2 after BuildNeighborList\n";
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LinkedCellCellManager2(args); }
void destroy(Plugin * m) { delete m; }

