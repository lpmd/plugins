//
//
//

#include <lpmd/physunits.h>

#include "ewald-2.h"

using namespace lpmd;

Ewald2::Ewald2(std::string args): Module("ewald2")
{
 ProcessArguments(args);
 status = true;
}

/*Funcion que asigna los valores por default principales de ewald
 */
void Ewald2::ReSet(SimulationCell &sc)
{
 kpoints.clear();
 double V = sc.Volume();
 long n = sc.Size();
 alpha = sqrt(M_PI)*powl(5.5*n/(V*V),0.166666667);
 SE = SelfEnergy(sc);
 rc = sqrt(11.5)/alpha;
 kc = 2*alpha*sqrt(11.5);
 Vector a1 = sc.Cell::GetVector(0);
 Vector a2 = sc.Cell::GetVector(1);
 Vector a3 = sc.Cell::GetVector(2);
 Vector b1 = (2*M_PI/(Dot(a1,Crux(a2,a3))))*Crux(a2,a3);
 Vector b2 = (2*M_PI/(Dot(a2,Crux(a3,a1))))*Crux(a3,a1);
 Vector b3 = (2*M_PI/(Dot(a3,Crux(a1,a2))))*Crux(a1,a2);
 int KMAX = 7;
 for(long kx = 0 ; kx <= KMAX ;kx++)
 {
  for(long ky = (kx==0 ? 0 : -KMAX)  ; ky <= KMAX ; ky++)
  {
   for(long kz = (kx==0 && ky==0 ? 1 : -KMAX); kz <= KMAX ; kz++)
   {
    Vector k = kx*b1 + ky*b2 + kz*b3;
    if(k.Mod()<kc) kpoints.push_back(k);
   }
  }
 }
 std::cerr << "************************************EWALD SETTINGS INFO**************************************" <<'\n';
 std::cerr << "Volumen = " << V << '\n';
 std::cerr << "N = " << n << '\n';
 std::cerr << "Vectores Base : " << '\n' << a1 << '\n' << a2 << '\n'<<a3 <<'\n';
 std::cerr << "Vectores Reciprocos : " << '\n' << b1 << '\n' << b2 << '\n' << b3 << '\n';
 std::cerr << "alpha = " << alpha << '\n';
 std::cerr << "ep = "<<ep <<'\n';
 std::cerr << "rc = " << rc << '\n';
 std::cerr << "kc = " << kc << '\n';
 std::cerr << "Numero de k-points en la suma = " << kpoints.size() << '\n';
 std::cerr << "Self-Energy  = " << SE << " [eV] " << '\n';
 std::cerr << "*********************************************************************************************" << '\n';
}

/* Calulo de las Energias, en primer lugar se implementa el metodo
 * virtual de Potential.
 */
double Ewald2::energy(SimulationCell &sc) 
{
 if(status==true) {ReSet(sc);status=false;}
 double RE = RealEnergy(sc);
 double KE = ReciprocalEnergy(sc);
 double DC = DipoleCorrectionEnergy(sc);
 std::cerr << "Real = " << RE << " Recip = " << KE << " Dipole = " << DC << " Self = " << SE << '\n';
 return RE+KE+DC+SE;
}
/*Suma la parte Real de la Energia, esta se suma como un sistema periodico (ver Rapaport SE).
 * Ya que no es necesario hacer la suma sobre las replicas, para eso se utiliza un rcut 
 * esferico en el sistema.
 * */
double Ewald2::RealEnergy(SimulationCell & sc) const
{
 double e=0.0e0;
 long n = sc.Size();
 for(int i=0;i<n;++i)
 {
  double qi = (sc.GetAtom(i)).Charge();
  for(int j=i+1;j<n;++j)
  {
   double qj = (sc.GetAtom(j)).Charge();
   double r = sc.Distance(i,j);
   if(r<rc)
   {
    e+=qi*qj*(1/r)*erfc(alpha*r);
   }
  }
 }
 return Q2a2EV*e; 
}
/*Suma la contribucion reciproca de la Energia
 */
double Ewald2::ReciprocalEnergy(SimulationCell & sc) const
{
 double e=0.0e0;
 long n = sc.Size();
 double V = sc.Volume();
 double qcos2 = 0.0e0;
 double qsin2 = 0.0e0;
 double f = 0.0e0;
 for(std::vector<Vector>::const_iterator it=kpoints.begin();it!=kpoints.end();++it)
 {
  Vector k = *it;
  double km = k.Mod();
  f = (1/(km*km)) * exp (-km*km/(4*alpha*alpha));
  qcos2=0.0e0;
  qsin2=0.0e0;
  for(long i=0;i<n;++i)
  {
   double qi = (sc.GetAtom(i)).Charge();
   qcos2 += qi*cos(Dot(k,(sc.GetAtom(i)).Position()));
   qsin2 += qi*sin(Dot(k,(sc.GetAtom(i)).Position()));
  }
  qcos2 = qcos2*qcos2;
  qsin2 = qsin2*qsin2;
  e+=(4*M_PI/V)*f*(qcos2+qsin2);
 }
 return Q2a2EV*e;
}
/*Asigna la contribucion de la auto-energia Sum_i qi**2 *alpha/sqrt(Pi)
 */
double Ewald2::SelfEnergy(SimulationCell & sc) const
{
 double e=0.0e0;
 long n=sc.Size();
 for(int i=0;i<n;++i)
 {
  double qi = (sc.GetAtom(i)).Charge();
  e+=qi*qi;
 }
 e=e*alpha/sqrt(M_PI);
 return -Q2a2EV*e;
}
/*Calcula la Correccion dipolar de la energia
 */
double Ewald2::DipoleCorrectionEnergy(SimulationCell & sc) const
{
 double e=0.0e0;
 long n = sc.Size();
 Vector tmp(0.0,0.0,0.0);
 for(int i=0;i<n;++i)
 {
  double qi = (sc.GetAtom(i)).Charge();
  Vector ri = (sc.GetAtom(i)).Position();
  tmp = tmp + qi*ri;
 }
 e = 2*M_PI*Dot(tmp,tmp)/((1+2*ep)*sc.Volume());
 return Q2a2EV*e;
}

/*Calculo de las fuerzas, en primer lugar se implementa metodo
 * virtual de Potential.
 * */
void Ewald2::UpdateForces(SimulationCell & sc)
{
 if(status==true) {ReSet(sc);status=false;}
 Vector ff, acci;
 long n = sc.Size();
 for (long i=0;i<n;++i)
 {
  const Atom & atom_i = sc.GetAtom(i);
  double mi = atom_i.Mass();
  acci = atom_i.Acceleration();
//  std::cerr << "acci = " << acci << '\n';
  Vector Real = RealForce(sc,i);
  Vector Reci = ReciprocalForce(sc,i);
  Vector Dipo = DipoleCorrectionForce(sc,i);
  ff = Real + Reci + Dipo;
  sc.SetAcceleration(i, acci + ff*(1/mi));
  if(i==(int)(n/2))
  {
   std::cerr << "mass = " << mi << '\n';
   std::cerr << "acci = " << acci << " acci.mod = "<<acci.Mod()<< '\n';
   std::cerr << "real = " << Real << " real.mod = "<<Real.Mod()<<'\n';
   std::cerr << "reci = " << Reci << " reci.mod = "<<Reci.Mod()<<'\n';
   std::cerr << "dipo = " << Dipo << " dipo.mod = "<<Dipo.Mod()<<'\n';
   std::cerr << "ff   = " << ff   << "   ff.mod = "<<ff.Mod()  <<'\n';
   Vector tmp = acci + ff*(1/mi);
   std::cerr << "newf = " << tmp  << " newf.mod = "<<tmp.Mod() <<'\n';
   std::cerr << "atom = " << (sc.GetAtom(i)).Acceleration() << '\n';
  }
 }
}

/*Calculo de contribucion real de fuerzas sobre el atomo i
 */
Vector Ewald2::RealForce(SimulationCell & sc, int i) const
{
 Vector fi(0.0,0.0,0.0);
 long n = sc.Size();
 for(long j=0;j<n;++j)
 {
  if(i!=j)
  {
   double qj = (sc.GetAtom(j)).Charge();
   Vector rij = sc.VectorDistance(i,j);
   double modrij = rij.Mod();
   if(modrij<rc)
   {
    double f = erfc(alpha*modrij)/modrij + (2*alpha/sqrt(M_PI))*exp(-alpha*alpha*modrij*modrij);
    fi = fi + (qj*f/(modrij*modrij))*rij;
   }
  }
 }
 return Q2a2FORCE*(sc.GetAtom(i)).Charge()*fi;
}

/*Calculo de contribucion reciproca de fuerzas sobre el atomo i
 */
Vector Ewald2::ReciprocalForce(SimulationCell & sc, int i) const
{
 Vector fi(0.0,0.0,0.0);
 long n = sc.Size();
 double V = sc.Volume();
 double qi = (sc.GetAtom(i)).Charge();
 for(std::vector<Vector>::const_iterator it=kpoints.begin();it!=kpoints.end();++it)
 {
  Vector k = *it;
  double km = k.Mod();
  Vector f = (qi/(km*km)) * exp (-km*km/(4*alpha*alpha)) * k;
  double qcos=0.0e0;
  double qsin=0.0e0;
  for(long j=0;j<n;++j)
  {
   double qj = (sc.GetAtom(j)).Charge();
   qcos += qj*cos(Dot(k,(sc.GetAtom(j)).Position()));
   qsin += qj*sin(Dot(k,(sc.GetAtom(j)).Position()));
  }
  qcos = sin(Dot(k,(sc.GetAtom(i)).Position()))*qcos;
  qsin = cos(Dot(k,(sc.GetAtom(i)).Position()))*qsin;
  fi = fi + (8*M_PI/V)*(qcos-qsin)*f;
 }
 return Q2a2FORCE*fi;
}

/*Calculo de contribucion dipolar de fuerzas sobre el atomo i
 */
Vector Ewald2::DipoleCorrectionForce(SimulationCell & sc, int i) const
{
 Vector fi(0.0,0.0,0.0);
 long n = sc.Size();
 for(long j=0;j<n;j++)
 {
  double qj = (sc.GetAtom(j)).Charge();
  Vector rj = (sc.GetAtom(j)).Position();
  fi = fi + qj*rj;
 }
 return -Q2a2FORCE*(4*M_PI/((1+2*ep)*sc.Volume()))*(sc.GetAtom(i)).Charge()*fi;
}


/*Asignacion y settings de los parametros iniciales para el calculo de ewald
 * solo valores y setups de rutina.
 */
void Ewald2::SetParameter(std::string name)
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

void Ewald2::Show() const
{
 Module::Show();
 if(alpha!=0) {std::cout << "   alpha = " << alpha << '\n';}
 else {std::cout << " alpha = will be calculated automatic during the evaluation." << '\n';}
 std::cout << "   ep    = " << ep << '\n';
}

std::string Ewald2::Keywords() const { return "alpha ep"; }

// Esto se inlcuye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) {return new Ewald2(args);}
void destroy(Module * m) { delete m; }
