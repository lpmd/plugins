//
//
//

#include "ewald.h"

#include <lpmd/session.h>
#include <lpmd/configuration.h>

using namespace lpmd;

Ewald::Ewald(std::string args): Plugin("ewald","2.0") 
{ 
 ParamList & params = (*this);
 DefineKeyword("kmax", "7");
 DefineKeyword("surfacedipole", "false");
 DefineKeyword("alpha");
 DefineKeyword("etol");
 ProcessArguments(args); 
 alpha = double(params["alpha"]);
 etol = double(params["etol"]);
 surfdip = bool(params["surfacedipole"]);
 kmax = int(params["kmax"]);
 rcut = sqrt(-log(etol))/alpha;
 kcut = 2.0*alpha*sqrt(-log(etol));
 std::cerr << "DEBUG EWALD -> rcut = " << rcut << '\n';
 kpoints = NULL;
}

Ewald::~Ewald() { delete kpoints; }

void Ewald::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      El modulo activa el potencial electrostatico calculado usando suma de    \n";
 std::cout << "      Ewald.                                                                   \n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << " Example                                                                       \n";
 std::cout << " Cargando el Modulo :                                                          \n";
 std::cout << " use ewald                                                                     \n";
 std::cout << "     alpha 0.37                                                                \n";
 std::cout << "     etol  0.0001                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " Llamando al modulo :                                                          \n";
 std::cout << " potential ewald                                                               \n\n";
 std::cout << "      De esta forma activa el potencial de ewald entre todas las especies      \n";
}

void Ewald::BuildKPointMesh(Configuration & conf)
{
 BasicCell & cell = conf.Cell();
 //BasicParticleSet & atoms = conf.GetAtoms();
 std::cerr << "DEBUG kcut = " << kcut << '\n';
 Vector R[3], K[3];  // real and reciprocal vectors
 for (int q=0;q<3;++q) R[q] = cell[q];
 K[0] = (2*M_PI/(Dot(R[0],Cross(R[1],R[2]))))*Cross(R[1],R[2]);
 K[1] = (2*M_PI/(Dot(R[1],Cross(R[2],R[0]))))*Cross(R[2],R[0]);
 K[2] = (2*M_PI/(Dot(R[2],Cross(R[0],R[1]))))*Cross(R[0],R[1]);
 kpoints = new std::vector<Vector>;
 int nrep = kmax;
 std::cerr << "DEBUG using nrep = " << nrep << '\n';
 for (int pp=0;pp<=nrep;++pp)
  for (int qq=(pp == 0 ? 0 : -nrep);qq<=nrep;++qq)
   for (int rr=(pp == 0 && qq == 0 ? 1 : -nrep);rr<=nrep;++rr)
   {
    Vector k = pp*K[0]+qq*K[1]+rr*K[2];
    if ((fabs(k.Module()) > 1.0E-10) && (k.Module() < kcut)) kpoints->push_back(k);
   } 
 std::cerr << "DEBUG number of k-points: " << kpoints->size() << '\n';
}

void Ewald::RealSpace(Configuration & conf, Vector * forces, double & e)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 int nrep = 0;
 double ep, e0;
 Vector ff;
 const double Q2a2EV = GlobalSession["q2a2ev"];
 const double Q2a2FORCE = GlobalSession["q2a2force"];

 e0 = 1.1/Q2a2EV;
 while (1)
 {
  //
  ep = 0.0;
  for (long int i=0;i<atoms.Size();++i) forces[i] = Vector(0.0, 0.0, 0.0);
  for (int pp=-nrep;pp<nrep;++pp)
   for (int qq=-nrep;qq<nrep;++qq)
    for (int rr=-nrep;rr<nrep;++rr)
    {
     Vector n = pp*cell[0]+qq*cell[1]+rr*cell[2];
     for (long int i=0;i<atoms.Size()-1;++i)
     {
      const double qi = atoms[i].Charge();
      for (long int j=i+1;j<atoms.Size();++j)
      {
       if ((((pp == qq) && (qq == rr)) && (rr == 0)) && (i == j)) continue; // FIXME nunca se cumple, con la suma como esta def.
       else
       {
        const double qj = atoms[j].Charge();
        Vector rij = atoms[j].Position()-atoms[i].Position();
        double rmod2 = (rij+n).SquareModule();
        if (rmod2 < rcut*rcut)
        {
	 double rmod = sqrt(rmod2);
         ep += qi*qj*erfc(alpha*rmod)/rmod;
         ff = (rij+n)*(1.0/rmod);
         ff = ff*qi*qj*(erfc(alpha*rmod)/rmod + 2.0*(alpha/sqrt(M_PI))*exp(-alpha*alpha*rmod*rmod));
         ff = ff*(-1.0/rmod);
         forces[i] = forces[i] + ff*(1.0/atoms[i].Mass());
         forces[j] = forces[j] - ff*(1.0/atoms[j].Mass());
        }
       }
      }
     }
    }
  if (fabs(ep - e0) < etol/Q2a2EV) break;
  e0 = ep; 
  nrep++;
 }
 for (long int i=0;i<atoms.Size();++i) atoms[i].Acceleration() =  atoms[i].Acceleration()+forces[i]*Q2a2FORCE;
 e = ep*Q2a2EV;
}

void Ewald::ReciprocalSpace(Configuration & conf, Vector * forces, double & e)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 const double Q2a2EV = GlobalSession["q2a2ev"];
 const double Q2a2FORCE = GlobalSession["q2a2force"];
 double ep = 0.0;
 for (long int i=0;i<atoms.Size();++i) forces[i] = Vector(0.0, 0.0, 0.0);
 if (kpoints == NULL) BuildKPointMesh(conf);
 for (unsigned int nk=0;nk<kpoints->size();++nk)
 {
  const Vector & k = (*kpoints)[nk];
  double kfac = 4.0*M_PI*(1.0/k.SquareModule())*(1.0/cell.Volume())*exp(-k.SquareModule()/(4.0*alpha*alpha));
  double sumqcos = 0.0, sumqsin = 0.0;
  for (long int i=0;i<atoms.Size();++i)
  {
   double dkr = Dot(k, atoms[i].Position());
   double qi = atoms[i].Charge();
   sumqcos += qi*cos(dkr);
   sumqsin += qi*sin(dkr);
  }  
  for (long int i=0;i<atoms.Size();++i)
  {
   const double qi = atoms[i].Charge();
   const Vector ri = atoms[i].Position();
   forces[i] = forces[i] + k*2.0*qi*kfac*(sin(Dot(k, ri))*sumqcos-cos(Dot(k, ri))*sumqsin);
  }
  ep += kfac*(pow(sumqcos, 2.0)+pow(sumqsin, 2.0));
 }
 for (long int i=0;i<atoms.Size();++i)
  atoms[i].Acceleration() = atoms[i].Acceleration()+forces[i]*(Q2a2FORCE/atoms[i].Mass()); 
 e = ep*Q2a2EV;
}

void Ewald::SurfaceDipole(Configuration & conf, Vector * forces, double & e)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 const double Q2a2EV = GlobalSession["q2a2ev"];
 const double Q2a2FORCE = GlobalSession["q2a2force"];
 e = 0.0;
 Vector v(0.0, 0.0, 0.0), sf(0.0, 0.0, 0.0);
 for (long int i=0;i<atoms.Size();++i)
 {
  forces[i][0] = 0.0e0;
  forces[i][1] = 0.0e0;
  forces[i][2] = 0.0e0;
  sf = sf + atoms[i].Charge()*atoms[i].Position();
 }
 for (long int i=0;i<atoms.Size();++i)
 {
  v = v + atoms[i].Charge()*atoms[i].Position();
  forces[i] = forces[i] - 4.0*M_PI/(6.0*cell.Volume())*atoms[i].Charge()*sf;
 }
 for (long int i=0;i<atoms.Size();++i)
  atoms[i].Acceleration() = atoms[i].Acceleration()+forces[i]*Q2a2FORCE; 
 e = 4.0*M_PI*v.SquareModule()/(6.0*cell.Volume())*Q2a2EV;
}

double Ewald::EnergyConstantCorrection(Configuration & conf)
{
 BasicParticleSet & atoms = conf.Atoms();
 BasicCell & cell = conf.Cell();
 const double Q2a2EV = GlobalSession["q2a2ev"];
 double e = 0.0;
 for (long int i=0;i<atoms.Size();++i) e += pow(atoms[i].Charge(), 2.0);
 e *= -alpha*(1.0/sqrt(M_PI));
 double totch = 0.0;
 for (long int i=0;i<atoms.Size();++i) totch += atoms[i].Charge();
 e -= pow(totch, 2.0)*M_PI/(4.0*cell.Volume()*alpha*alpha);
 return e*Q2a2EV;
}

double Ewald::energy(Configuration & conf) 
{ 
 if (conf.HaveAny(Tag("pe"))) return double(Parameter(conf.GetTag(conf, Tag("pe"))));
 throw PluginError("ewald", "This shouldn\'t happen");
}

void Ewald::UpdateForces(Configuration & conf) 
{ 
 BasicParticleSet & atoms = conf.Atoms();
 double ereal=0.0, erecip=0.0, edip=0.0, ecorr=0.0;
 Vector * forces = new Vector[atoms.Size()];
 RealSpace(conf, forces, ereal);
 ReciprocalSpace(conf, forces, erecip); 
 if (surfdip) SurfaceDipole(conf, forces, edip);
 ecorr = EnergyConstantCorrection(conf); 
 /*
 std::cerr << "DEBUG Energy real part = " << ereal << '\n';
 std::cerr << "DEBUG Energy recip part = " << erecip << '\n';
 std::cerr << "DEBUG Energy surf dipole = " << edip << '\n';
 std::cerr << "DEBUG Energy const correct = " << ecorr << '\n';
 */
 double e = ereal+erecip+edip+ecorr;
 conf.SetTag(conf, Tag("pe"), e);
 delete [] forces;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Plugin * create(std::string args) { return new Ewald(args); }
void destroy(Plugin * m) { delete m; }


