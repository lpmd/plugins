//Finnis-Sinclair Extended Potential
//Dai, X. D., Kong, Y., Li, J. H., and Liu, B. X., 2006, J. Phys.: Condens. Matter, 18, 45274542. 30
//

#include "finnissinclair-ext.h"
#include <iostream>

using namespace lpmd;

FinnisSinclairExt::FinnisSinclairExt(std::string args): Plugin("finnissinclair-ext", "2.1")
{
 ParamList & params = (*this);
 ProcessArguments(args); 
 c0 = params["c0"];
 c1 = params["c1"];
 c2 = params["c2"];
 c3 = params["c3"];
 c4 = params["c4"];
 A = params["A"];
 B = params["B"];
 c = params["c"];
 d = params["d"];
 if(d>=c) SetCutoff(d);
 else SetCutoff(c);
}

void FinnisSinclairExt::ShowHelp() const
{
 std::cout << " General Info      >>                                                          \n";
 std::cout << "      This plugin incorporate the FinnisSinclairExt potential used frequently  \n";
 std::cout << " for metallic atomic interaction. Based in embedded atom model.                \n\n";
 std::cout << " V(r) = (r-c)^2*(c0+c1*r+c2*r^2+c3*r^3+c4*r^4) ; rho(r) = (r-d)^2+B^2*(r-d)^4 ;\n";
 std::cout << " F(rho) = -A/sqrt(rho)                                                         \n\n";
 std::cout << " General Options   >>                                                          \n";
 std::cout << "      c0            : Value of c0 for the potential.  [eV/A^2]                 \n";
 std::cout << "      c1            : Value of c1 for the potential.  [eV/A^3]                 \n";
 std::cout << "      c2            : Value of c2 for the potential.  [eV/A^4]                 \n";
 std::cout << "      c3            : Value of c3 for the potential.  [eV/A^5]                 \n";
 std::cout << "      c4            : Value of c4 for the potential.  [eV/A^6]                 \n";
 std::cout << "      A             : Value of  A for the potential.  [eV/A]                   \n";
 std::cout << "      B             : Value of  B for the potential.  [1/A^2]                  \n";
 std::cout << "      c             : First cutoff of the system (cutoff_1=c). [A]             \n";
 std::cout << "      d             : Second cutoff od the system (cutoff_2=d).[A]             \n";
 std::cout << " Example          >>                                                           \n";
 std::cout << " #Loading the plugin                                                           \n";
 std::cout << " use FinnisSinclairExt as FS                                                   \n";
 std::cout << "     c0  47.98066                                                               \n";
 std::cout << "     c1 −34.09924                                                               \n";
 std::cout << "     c2  5.832293                                                               \n";
 std::cout << "     c3  0.017494                                                               \n";
 std::cout << "     c4  0.020393                                                               \n";
 std::cout << "     A   1.848648                                                               \n";
 std::cout << "     B   0.000000                                                               \n";
 std::cout << "     c   3.257200                                                               \n";
 std::cout << "     d   4.147200                                                              \n";
 std::cout << " enduse                                                                        \n";
 std::cout << " #using the loaded plugin                                                      \n";
 std::cout << " potential FS Mo Mo                                                            \n";
 std::cout << " #Values from 2006, J. Phys.: Condens. Matter, 18, 45274542. 30              \n\n";
}

double FinnisSinclairExt::pairEnergy(const double &r) const
{
	if(r<=c) 
    {
     double r2 = r*r;
     return (r-c)*(r-c)*(c0+c1*r+c2*r2+c3*r2*r+c4*r2*r2);
    }
    else return 0.0e0;
}

double FinnisSinclairExt::rhoij(const double &r) const
{
	if(r<=d)
    {
     double rd=r-d; 
     return (rd*rd)*(1.0e0+B*B*rd*rd);
    }
    else return 0.0e0;
}

double FinnisSinclairExt::F(const double &rhoi) const
{
	return -A*sqrt(rhoi);
}

Vector FinnisSinclairExt::PairForce(const Vector &normrij, const double &r) const
{
    if(r<=c)
    {
     double r2=r*r;
     double rc=r-c;
	 double t1=2.0e0*rc*(c0+c1*r+c2*r2+c3*r2*r+c4*r2*r2);
	 double t2=(c1+2.0e0*c2*r+3.0e0*c3*r2+4.0e0*c4*r2*r)*rc*rc;
	 return (t1+t2)*normrij;
    }
    else return Vector(0.0e0,0.0e0,0.0e0);
}

Vector FinnisSinclairExt::ManyBodies(const Vector &normrij, const double &rhoi, const double &rhoj, const double &r) const
{
    if(r<=d)
    {
     double rd=r-d;
	 double t1=A*0.5e0*((1/sqrt(rhoi))+(1/sqrt(rhoj)));
	 double t2=2.0e0*rd + 4.0e0*B*B*rd*rd*rd;
	 return -t1*t2*normrij;
    }
    else return Vector(0.0e0,0.0e0,0.0e0);
}

// Esto se incluye para que el modulo pueda ser cargado dinámicamente
Plugin * create(std::string args) {return new FinnisSinclairExt(args);}
void destroy(Plugin * m) { delete m; }

