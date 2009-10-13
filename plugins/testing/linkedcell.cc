//
//
//

#include "linkedcell.h"

using namespace lpmd;

inline int WrapAround(int i, int n)
{
 if (i < 0) i += ((-i/n)+1)*n;
 else if (i > n-1) i -= (i/n)*n;
 return i;
}

LinkedCell::LinkedCell(std::string args): Plugin("linkedcell", "1.0")
{ 
 ParamList & params = (*this);
 DefineKeyword("cutoff", "7.0");
 DefineKeyword("nx", "0");
 DefineKeyword("ny", "0");
 DefineKeyword("nz", "0");
 DefineKeyword("mode", "noauto");
 // 
 ProcessArguments(args);
 cutoff = double(params["cutoff"]);
 nx = int(params["nx"]);
 ny = int(params["ny"]);
 nz = int(params["nz"]);
 if (nx==0 || ny==0 || nz==0) mode=true;
 else if((nx>=4 && ny>=4) && nz>=4) mode=false;
 else mode=true;
 //
 if (params["mode"]=="auto" || params["mode"]=="AUTO") mode=true;
 else mode=false;
 // 
 //
 //
 head = tail = 0;
 atomlist = indexc = 0;
 subcell = 0;
}

LinkedCell::~LinkedCell() 
{ 
 delete [] head;
 delete [] tail;
 delete [] atomlist;
 delete [] subcell;
 delete [] indexc;
}

void LinkedCell::Reset() { }

void LinkedCell::UpdateCell(Configuration & conf) 
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 if (atomlist == 0) atomlist = new long[atoms.Size()];
 if (indexc == 0) indexc = new long[atoms.Size()];
 if (mode == true)
 {
  double minx = cell[0].Module();
  double miny = cell[1].Module();
  double minz = cell[2].Module();
  double fnn = (cutoff/2.5); //NOTE : Approximated value.
  int n = 1;
  while(true)
  {
   double previo = fabs(minx/n - fnn);
   n++;
   double actual = fabs(minx/n - fnn);
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
        z[0] = WrapAround(i+dx, nx);
        z[1] = WrapAround(j+dy, ny);
        z[2] = WrapAround(k+dz, nz);
        long p = z[2]*(nx*ny)+z[1]*nx+z[0];
        subcell[q*cells_inside+(r++)] = p;
       }
    }
 }

 //
 for (long q=0;q<nx*ny*nz;++q) head[q] = tail[q] = -1;
 for (long r=0;r<atoms.Size();++r)
 {
  const Vector fpos = cell.Fractional(atoms[r].Position());
  //
  int i = WrapAround(int(floor(nx*fpos[0])), nx);
  int j = WrapAround(int(floor(ny*fpos[1])), ny);
  int k = WrapAround(int(floor(nz*fpos[2])), nz);
  long q = k*(nx*ny)+j*nx+i;
  indexc[r] = q;
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
// const Vector fpos = cell.Fractional(atoms[i].Position());
 //
// int p = WrapAround(int(floor(nx*fpos[0])), nx);
// int q = WrapAround(int(floor(ny*fpos[1])), ny);
// int r = WrapAround(int(floor(nz*fpos[2])), nz);
 long cind = indexc[i];//r*(nx*ny)+q*nx+p;
 //
 AtomPair nn;
 nlist.Clear();
 nn.i = &atoms[i];
 nn.i_index = i;
 for (int c=0;c<cells_inside;++c)
 {
  int neighbor_cell = subcell[cind*cells_inside+c];
  for (long z=head[neighbor_cell];z != -1;z=atomlist[z])
  {
   if (z == i) continue;
   if ((full == false) && (z > i)) continue;
   nn.j = &atoms[z];
   nn.rij = cell.Displacement(nn.i->Position(), nn.j->Position());
   nn.r2 = nn.rij.SquareModule();
   nn.j_index = z;
   if (nn.r2 < rcut*rcut) nlist.Append(nn);
  }
 }
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new LinkedCell(args); }
void destroy(Plugin * m) { delete m; }

