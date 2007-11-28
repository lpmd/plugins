//
//
//

#include <lpmd/physunits.h>

#include "ewald.h"

using namespace lpmd;

Ewald::Ewald(std::string args): Module("ewald")
{
 ProcessArguments(args); 
}

double Ewald::RealEnergy(SimulationCell & sc) const
{
 double e = 0.0;
 double rc = sqrt(11.5)/alpha;
 long n = sc.Size();
 for (long nx=-1;nx<=1;nx++)  
 {
  for (long ny=-1;ny<=1;ny++)
  {
   for (long nz=-1;nz<=1;nz++)
   {
    for (long i=0;i<=n;++i)
    {
     for (long j=i+1;j<n;++j)
     {
      double qi = (sc.GetAtom(i)).Charge();
      double qj = (sc.GetAtom(j)).Charge();
      if(nx==0 && (ny == 0 && nz ==0))
      {
       double r=sc.RealDistance(i,j);
       if(i !=j && r<rc) e+=qi*qj*erfc(alpha*r)/r;
      }
      else
      {
       Vector n=nx*sc.Cell::GetVector(0) + ny*sc.Cell::GetVector(1) + nz*sc.Cell::GetVector(2);
       Vector rij = sc.VectorRealDistance(i,j);
       double r = (rij+n).Mod();
       if(r<rc) e+=qi*qj*erfc(alpha*r)/r;
      }
     }
    }
   }
  }
 }  		 
 return e*Q2a2EV;
}

double Ewald::ReciprocalEnergy(SimulationCell &sc) const
{
 double e = 0.0e0;
 double kc = 2*alpha*sqrt(11.5);
 //std::cerr << "kc = " << kc << '\n';
 long n = sc.Size();
 double V = sc.Volume();
 int counter=0;
 for(long kx=0;kx<=7;kx++)
 {
  for(long ky=-7;ky<=7;ky++)
  {
   for(long kz=-7;kz<=7;kz++)
   {
    Vector a1 = sc.Cell::GetVector(0);
    Vector a2 = sc.Cell::GetVector(1);
    Vector a3 = sc.Cell::GetVector(2);
    Vector k1 = (2*M_PI/(Dot(a1,Crux(a2,a3))))*kx*Crux(a2,a3);
    Vector k2 = (2*M_PI/(Dot(a2,Crux(a3,a1))))*ky*Crux(a3,a1);
    Vector k3 = (2*M_PI/(Dot(a3,Crux(a1,a2))))*kz*Crux(a1,a2);
    Vector k = k1+k2+k3;
//    Vector k(2*M_PI*kx/V,2*M_PI*ky/V,2*M_PI*kz/V);
    if(k.Mod()<kc && (ky!=0 && kz!=0))
    {
     double tmp=(1/k.Mod2())*exp(-k.Mod2()/(4*alpha*alpha));
     double cos1=0.0e0;
     double sin1=0.0e0;
     for(long j=0;j<n;++j)
     {
      double qj = (sc.GetAtom(j)).Charge();
      cos1 += qj*cos(Dot(k,(sc.GetAtom(j)).Position()));
      sin1 += qj*sin(Dot(k,(sc.GetAtom(j)).Position()));
     }
     cos1 = fabs(cos1)*fabs(cos1);
     sin1 = fabs(sin1)*fabs(sin1);
     e+=tmp*(cos1+sin1);
     counter++;
    }
   }
  }
 }
 std::cerr << " Evaluated with " << counter << " k-points." << '\n';
 return e*Q2a2EV*(4*M_PI/V);
}

double Ewald::SelfEnergy(SimulationCell &sc) const
{
 double e=0.0e0;
 long n = sc.Size();
 for(int i=0;i<n;++i)
 {
  double a=(sc.GetAtom(i)).Charge();
  e += a*a ;
 }
 e = -(alpha/sqrt(M_PI))*e;
 return e*Q2a2EV;
}

double Ewald::DipoleCorrectionEnergy(SimulationCell &sc) const
{
 double e=0.0e0;
 long n = sc.Size();
 double V = sc.Volume();
 Vector tmp(0,0,0);
 for(long i=0;i<n;++i)
 {
  double qi = (sc.GetAtom(i)).Charge();
  tmp = tmp + qi * (sc.GetAtom(i)).Position();
 }
 e = tmp.Mod2();
 return Q2a2EV*(2*M_PI/((1+2*ep)*V))*e;
}

double Ewald::energy(SimulationCell &sc) 
{
 double RE = RealEnergy(sc);
 double ReE= ReciprocalEnergy(sc);
 double SE = SelfEnergy(sc);
 double DE = DipoleCorrectionEnergy(sc);
 std::cerr << "REAL = " << RE << " \t RECIPROCAL = " << ReE << " \t SELF = " << SE << " \t DIPOLE = " << DE << '\n';
 return RE + ReE + SE + DE;
 //return RealEnergy(sc) + ReciprocalEnergy(sc) + SelfEnergy(sc) + DipoleCorrectionEnergy(sc);
}

Vector Ewald::RealForce(SimulationCell &sc, int i) const
{
 long n = sc.Size();
 double rc = sqrt(11.5)/alpha;
 Vector Fi(0,0,0);
 double qi=(sc.GetAtom(i)).Charge();
 for(long nx = -1;nx<=1;nx++)
 {
  for(long ny = -1;ny<=1;ny++)
  {
   for(long nz = -1;nz<=1;nz++)
   {
    for(long j = 0;j<n;++j)
    {
     if(nx==0 && (ny ==0 && nz == 0))
     {
      if(j!=i)
      {
       double qj = (sc.GetAtom(j)).Charge();
       double r = sc.RealDistance(i,j);
       if(r<rc)
       {
	double t = (2*alpha/sqrt(M_PI))*exp(-alpha*alpha*r*r) + erfc(alpha*r)/r;
	Vector Fij = sc.VectorRealDistance(i,j);
	Fij.Scale(qi*qj*t/(r*r));
	Fi=Fi+Fij;
       }
      }
     }
     else
     {
      Vector rij = sc.VectorRealDistance(i,j);
      Vector n = nx*sc.Cell::GetVector(0) + ny*sc.Cell::GetVector(1) + nz*sc.Cell::GetVector(2);
      double qj = (sc.GetAtom(j)).Charge();
      Vector rijpn = rij + n ;
      double r = rijpn.Mod();
      if(r<rc)
      {
       double t = (2*alpha/sqrt(M_PI))*exp(-alpha*alpha*r*r) + erfc(alpha*r)/r;
       //Vector Fij = sc.VectorRealDistance(i,j);
       rijpn.Scale(qi*qj*t/(r*r));
       Fi=Fi+rijpn;
      }
     }
    }
   }
  }
 }
 return Q2a2FORCE*Fi;
}

Vector Ewald::ReciprocalForce(SimulationCell &sc, int i) const
{
 long n = sc.Size();
 double kc = 2*alpha*sqrt(11.5);
 double V = sc.Volume();
 Vector Fi(0,0,0);
 double ff = 0.0e0;
 double qi = (sc.GetAtom(i)).Charge();
 for(long kx=0;kx<=7;kx++)
 {
  for(long ky=-7;ky<=7;ky++)
  {
   for(long kz=-7;kz<=7;kz++)
   {
    Vector a1 = sc.Cell::GetVector(0);
    Vector a2 = sc.Cell::GetVector(1);
    Vector a3 = sc.Cell::GetVector(2);
    Vector k1 = (2*M_PI/(Dot(a1,Crux(a2,a3))))*kx*Crux(a2,a3);
    Vector k2 = (2*M_PI/(Dot(a2,Crux(a3,a1))))*ky*Crux(a3,a1);
    Vector k3 = (2*M_PI/(Dot(a3,Crux(a1,a2))))*kz*Crux(a1,a2);
    Vector k = k1+k2+k3;
//    Vector k(2*M_PI*kx/V,2*M_PI*ky/V,2*M_PI*kz/V);
    if(k.Mod()<kc && (ky!=0 && kz!=0))
    {
     Vector ri = (sc.GetAtom(i)).Position();
     double k2 = k.Mod2();
     double factor = qi*(1/k2)*exp(-k2/(4*alpha*alpha));
     double sum1=0.0e0;
     double sum2=0.0e0;
     for(long j=0;j<n;++j)
     {
      double qj = (sc.GetAtom(j)).Charge();
      Vector rj = (sc.GetAtom(j)).Position();
      sum1+=qj*cos(Dot(k,rj));
      sum2+=qj*sin(Dot(k,rj));
     }
     sum1 = sin(Dot(k,ri))*sum1;
     sum2 = cos(Dot(k,ri))*sum2;
     
     ff+=factor*(sum1 - sum2);
     k.Scale(ff);
     Fi = Fi + k;
    }
   }
  }
 }
 return Q2a2FORCE*Fi*(8*M_PI/V);
}

Vector Ewald::DipoleCorrectionForce(SimulationCell &sc, int i) const
{
 long n = sc.Size();
 double V = sc.Volume();
 Vector Fi(0,0,0);
 //return Fi;
 double qi = (sc.GetAtom(i)).Charge();
 for(long j = 0; j<n ; ++j)
 {
  double qj = (sc.GetAtom(j)).Charge();
  Fi = Fi + qi*qj*(sc.GetAtom(j)).Position();
 }
 return Q2a2FORCE*(2*M_PI/((1+2*ep)*V))*Fi;
}

void Ewald::UpdateForces(SimulationCell & sc)
{ 
 Vector ff, acci;
 long n = sc.Size();
 Vector RealAv(0,0,0);
 Vector ReciAv(0,0,0);
 Vector DipoAv(0,0,0);
 for (long i=0;i<n;++i)
 {
  double mi = sc.GetAtom(i).Mass();
  acci = sc.GetAtom(i).Acceleration();
  Vector Real = RealForce(sc,i);
  Vector Reci = ReciprocalForce(sc,i);
  Vector Dipo = DipoleCorrectionForce(sc,i);
  ff = Real + Reci + Dipo;
  RealAv = RealAv + Real;
  ReciAv = ReciAv + Reci;
  DipoAv = DipoAv + Dipo;
//  std::cerr << "RealForce = " << RealForce(sc,i) << "\t ReciprocalForce = " << ReciprocalForce(sc,i) << "\t DipoleForce = " << DipoleCorrectionForce(sc,i) << '\n';
//  std::cerr << "Total Force = " << ff << '\n';
//  std::cerr << "Previous acci = " << acci << "\t New Accel = " << acci + ff*(1/mi) << '\n';
  sc.SetAcceleration(i, acci + ff*(FORCEFACTOR/mi));
 }
// std::cerr << " RealForce Average = " << RealAv << " \t and RealAv/n = " << RealAv/n << '\n';
// std::cerr << " ReciForce Average = " << ReciAv << " \t and ReciAv/n = " << ReciAv/n << '\n';
// std::cerr << " DipoForce Average = " << DipoAv << " \t and DipoAv/n = " << DipoAv/n << '\n';
}

void Ewald::SetParameter(std::string name)
{
 if(name == "alpha")
 {
  AssignParameter("alpha",GetNextWord());
  alpha = GetDouble("alpha");
 }
 if(name == "ep")
 {
  AssignParameter("ep",GetNextWord());
  ep = GetDouble("ep");
 }
}

void Ewald::Show() const
{
 Module::Show();
 std::cout << "   alpha = " << alpha << '\n';
 std::cout << "   ep    = " << ep << '\n';
}

std::string Ewald::Keywords() const { return "alpha ep"; }

// Esto se inlcuye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) {return new Ewald(args);}
void destroy(Module * m) { delete m; }
