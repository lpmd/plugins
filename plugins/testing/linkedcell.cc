//
//
//

#include "linkedcell.h"

using namespace lpmd;

LinkedCell::LinkedCell(std::string args): Plugin("linkedcell", "2.0")
{ 
 ParamList & params = (*this);
 DefineKeyword("cutoff", "7.0");
 DefineKeyword("nx", "6");
 DefineKeyword("ny", "6");
 DefineKeyword("nz", "6");
 DefineKeyword("mode", "noauto");
 // 
 ProcessArguments(args);
 cutoff = double(params["cutoff"]);
 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);
 if ((nx==0 || ny==0) || nz==0) mode=true;
 else if((nx>=1 && ny>=1) && nz>=1) mode=false;
 if (params["mode"]=="auto" || params["mode"]=="AUTO") mode=true;
 // 
 //
 //
 head = tail = 0;
 atomlist = 0;
 subcell = 0;
}

LinkedCell::~LinkedCell() 
{ 
 delete [] head;
 delete [] tail;
 delete [] atomlist;
 delete [] subcell;
}

void LinkedCell::Reset() { }

void LinkedCell::UpdateCell(Configuration & conf) 
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 if (atomlist == 0) atomlist = new long[atoms.Size()];
 if (mode == true)
 {
  double minx = cell[0].Module();
  double miny = cell[1].Module();
  double minz = cell[2].Module();
  double fnn = double(cutoff/5.0e0); //NOTE : Approximated value. one atom by cell.
  int n = 1;
  while(true)
  {
   double previo = fabs(minx/n - fnn);
   n++;
   double actual = fabs(minx/n-fnn);
   if (actual<1E-3) break;
   else if (actual>previo) {n--;break;}
   else continue;
  }
  if ((n%2)!=0) n--;
  nx = n;
  DebugStream() << "-> Using nx = " << nx << " subdivision scheme."<<'\n';
  n=1;
  while(true)
  {
   double previo = fabs(miny/n - fnn);
   n++;
   double actual = fabs(miny/n - fnn);
   if(actual<1E-3) break;
   else if (actual>previo) {n--;break;}
   else continue;
  }
  if ((n%2)!=0) n--;
  ny = n;
  DebugStream() << "-> Using ny = " << ny << " subidvision scheme."<<'\n';
  n=1;
  while(true)
  {
   double previo = fabs(minz/n - fnn);
   n++;
   double actual = fabs(minz/n - fnn);
   if(actual<1E-3) break;
   else if (actual>previo) {n--;break;}
   else continue;
  }
  if ((n%2)!=0) n--;
  nz = n;
  DebugStream() << "-> Using nz = " << nz << " subdivision scheme." << '\n';
  n=1;
  mode=false;
 }
 //
 //
 //

 if (head != 0) delete [] head;
 head = new int[nx*ny*nz];
 if (tail != 0) delete [] tail;
 tail = new int[nx*ny*nz];
 //
 if (subcell == 0)
 {
  double d = cell[0].Module()/double(nx);
  if (cell[1].Module()/double(ny) < d) d = cell[1].Module()/double(ny);
  if (cell[2].Module()/double(nz) < d) d = cell[2].Module()/double(nz);
  int side = int(ceil((cutoff/d)-0.5));
  cells_inside = (2*side+1)*(2*side+1)*(2*side+1);
  if (cells_inside < 27 || d > cutoff) 
  {
   nx+=2;ny+=2;nz+=2;
   DebugStream() << "-> Extremely small values of nx,ny,nz or d>cutoff."<<'\n';
   DebugStream() << "-> Update nx-ny-nz to = " << nx <<"-"<<ny<<"-"<<nz<<'\n';
   if (subcell != 0) { delete [] subcell; subcell = 0; }
   if (atomlist !=0 ) { delete [] atomlist; atomlist = 0; }
   //nwin=nfail=0;
   UpdateCell(conf);
  }
  DebugStream() << "-> Using " << cells_inside << " neighboring subcells\n";
  subcell = new int[(nx*ny*nz)*cells_inside];
  int z[3];
  for (int k=0;k<nz;++k)
   for (int j=0;j<ny;++j)
    for (int i=0;i<nx;++i)
    {
     int r = 0;
     long q = k*(nx*ny)+j*nx+i;
     for (int dz=-side;dz<=side;++dz)
      for (int dy=-side;dy<=side;++dy)
       for (int dx=-side;dx<=side;++dx)
       {
        z[0] = i+dx;
        z[1] = j+dy;
        z[2] = k+dz;
        if (z[0] < 0) z[0] += nx;
        else if (z[0] > nx-1) z[0] -= nx;
        if (z[1] < 0) z[1] += ny;
        else if (z[1] > ny-1) z[1] -= ny;
        if (z[2] < 0) z[2] += nz;
        else if (z[2] > nz-1) z[2] -= nz;
        long p = z[2]*(nx*ny)+z[1]*nx+z[0];
        subcell[q*cells_inside+(r++)] = p;
       }
    }
 }

 //
 for (long q=0;q<nx*ny*nz;++q) head[q] = tail[q] = -1;
 long r = 0;
#ifdef _OPENMP
#pragma omp parallel for private ( r )
#endif
 for (r=0;r<atoms.Size();++r)
 {
  const Vector fpos = cell.Fractional(atoms[r].Position());
  //
  int i = int(floor(nx*fpos[0]));
  int j = int(floor(ny*fpos[1]));
  int k = int(floor(nz*fpos[2]));
//  if (i < 0) i += nx;
  if (i > nx-1) i -= nx;
//  if (j < 0) j += ny;
  if (j > ny-1) j -= ny;
//  if (k < 0) k += nz;
  if (k > nz-1) k -= nz;
  long q = k*(nx*ny)+j*nx+i;
  //
  if (head[q] == -1) head[q] = tail[q] = r;
  else 
  {
   atomlist[tail[q]] = r;
   tail[q] = r; 
  }
  atomlist[r] = -1;
 }
}

void LinkedCell::BuildNeighborList(Configuration & conf, long i, NeighborList & nlist, bool full, double rcut)
{ 
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 const Vector fpos = cell.Fractional(atoms[i].Position());
 //
 int p = int(floor(nx*fpos[0]));
 int q = int(floor(ny*fpos[1]));
 int r = int(floor(nz*fpos[2]));
// if (p < 0) p += nx;
 if (p > nx-1) p -= nx;
// if (q < 0) q += ny;
 if (q > ny-1) q -= ny;
// if (r < 0) r += nz;
 if (r > nz-1) r -= nz;
 long cind = r*(nx*ny)+q*nx+p;
 //
 AtomPair nn;
 nlist.Clear();
 nn.i = &atoms[i];
 int * c=0;long z=0;
 c = &(subcell[cind*cells_inside]);
#ifdef _OPENMP
#pragma omp parallel for
#endif
 for (int q=0;q<cells_inside;++q)
 {
  z = atomlist[*(c++)];
  if ((z < 0) || (z == i)) continue;
  if ((z > i) && (full == false)) continue;
  nn.j = &atoms[z];
  nn.rij = cell.Displacement(nn.i->Position(), nn.j->Position());
  nn.r2 = nn.rij.SquareModule();
  if (nn.r2 < rcut*rcut) nlist.Append(nn);
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LinkedCell(args); }
void destroy(Plugin * m) { delete m; }

