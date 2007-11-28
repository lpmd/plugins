//
//
//


#include <sstream>

#include <lpmd/util.h>

#include "povray.h"

using namespace lpmd;

class POVCylinder
{
 public:
  POVCylinder(Vector a,Vector b,double r,Vector c){base=a;cap=b;rad=r;col=c;}
  ~POVCylinder(){};
  std::string GetCode()
  {
   std::ostringstream ostr;
   ostr << "cylinder {\n";
   ostr << base.GetPOV() <<", "<<cap.GetPOV()<<", "<<rad<<"\n";
   ostr << "         texture { \n";
   ostr << "                    pigment {color "<<col.GetPOV()<<"}\n";
   ostr << "                 }\n";
   ostr << "}\n\n";
   return ostr.str();
  }
 private:
  Vector base;
  Vector cap;
  Vector col;
  double rad;
};

class POVSphere
{
 public:
  POVSphere(Vector a,double r,Vector c){pos=a;rad=r;col=c;}
  ~POVSphere(){};
  std::string GetCode()
  {
   std::ostringstream ostr;
   ostr << "sphere { \n";
   ostr << "\t " << pos.GetPOV() << ", "<<rad*0.75<<"\n";
   ostr << "   texture  { \n";
   ostr << "                pigment {color "<<col.GetPOV()<<"}\n";
   ostr << "            } \n";
   ostr << "}" << '\n';
   return ostr.str();
  }
 private:
  Vector pos,col;
  double rad;
};

class POVText
{
 public:
  POVText(std::string t,Vector p,Vector r,double x,double y, Vector d){text=t;pos=p;rot=r;a=x,b=y;col=d;}
  ~POVText(){};
  std::string GetCode()
  {
   std::ostringstream ostr;
   ostr << "text { \n";
   ostr << "\t ttf\"timrom.ttf\" \"" << text << "\" "<<a<<", "<<b<<"\n";
   ostr << "\t pigment {color "<<col.GetPOV()<<"}\n";
   ostr << "\t rotate    "<<rot.GetPOV()<<"\n";
   ostr << "\t translate "<<pos.GetPOV()<<"\n";
   ostr << "}" << '\n';
   return ostr.str();
  }
 private:
  Vector pos,col,rot;
  double a,b;
  std::string text;
};

POVRAY::POVRAY(std::string args): Module("povray")
{
 direct = "povray";
 header = "movie-";
 counter = 0;
 level=0;
 ProcessArguments(args);
}

POVRAY::~POVRAY()
{
}

void POVRAY::SetParameter(std::string name)
{
 if (name == "start")
 {
  AssignParameter("start", GetNextWord());
  start_step = GetInteger("start");
 }
 if (name == "end")
 {
  AssignParameter("end", GetNextWord());
  end_step = GetInteger("end");
 }
 if (name == "each")
 {
  AssignParameter("each", GetNextWord());
  interval = GetInteger("each");
 }
 if (name == "header")
 {
  AssignParameter("header", GetNextWord());
  header = GetString("header");
 }
 if (name == "direct")
 {
  AssignParameter("direct", GetNextWord());
  direct = GetString("direct");
 }
}

void POVRAY::Show() const
{
 Module::Show();
 std::cout << "   start    = " << start_step << '\n';
 std::cout << "   end      = " << end_step << '\n';
 std::cout << "   step     = " << interval << '\n';
 std::cout << "   header   = " << header << '\n';
 std::cout << "   direct   = " << direct << '\n';
}

void POVRAY::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = POVRAY2 " << '\n';
 std::cout << " Module Version     = 1.0 " << '\n';
 std::cout << " Support API lpmd   = 1.0 " << '\n';
 std::cout << " Problems Report to = gnm@gnm.cl " << '\n';
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string POVRAY::Keywords() const { return "start end step header direct"; }

void POVRAY::Apply(const MD & md)
{
 SimulationCell sc=md.GetCell();
 std::ostringstream ostr;
 Vector a = sc.GetVector(0);
 Vector b = sc.GetVector(1);
 Vector c = sc.GetVector(2);
 ostr << "mkdir -p " << direct << ";";
 system (ostr.str().c_str());
 ostr.str("");
 ostr <<direct<<"/"<< header ;
 ostr.setf(std::ios::right,std::ios::adjustfield);
 ostr.fill('0');
 ostr << std::setw(10) << counter << ".pov";
 std::ofstream output (ostr.str().c_str(),std::ios::out);
 //Set basic file headers.
 output << "#include \"colors.inc\"" << '\n';
 output << "#include \"stones.inc\"" << '\n';
 output << "#include \"textures.inc\"" << '\n';
 output << "#include \"shapes.inc\"" << '\n';
 output << "#include \"glass.inc\"" << '\n';
 output << "#include \"metals.inc\"" << '\n';
 output << "#include \"woods.inc\"" << '\n';
 output << "global_settings { ambient_light rgb<1, 1, 1> }" << '\n';
 //Set camera and lights.
 Vector camerapos = 1.25*a + 1.25*b -2*c;
 Vector cameralok = (a+b+c)/2;
 Vector direction = camerapos-cameralok;
 double angx = Ang(direction,a);
 double angy = Ang(direction,b);
 Vector corner1(0,0,0);
 Vector corner2 = a+b+c;
 Vector light1(camerapos.GetX(),camerapos.GetY(),camerapos.GetZ());
 output << "camera {\n";
 output << "\t location " << camerapos.GetPOV() << '\n';
 output << "\t look_at  " << cameralok.GetPOV() << '\n';
 output << "}\n\n"; 
 output << "light_source { " << light1.GetPOV() << "color White}"<< '\n';
 //Set Cilindros
 Vector p0(0,0,0);
 Vector p1=a,p2=b,p3=c;
 Vector p4=p1+p2,p5=p1+p3,p6=p3+p2,p7=p1+p2+p3;
 Vector ccol(1,0,0);
 double rad = ElemRad[sc.GetAtom(0).Species()];
 POVCylinder C1(p0,p1,rad*0.1,ccol);output << C1.GetCode();
 POVCylinder C2(p0,p2,rad*0.1,ccol);output << C2.GetCode();
 POVCylinder C3(p0,p3,rad*0.1,ccol);output << C3.GetCode();
 POVCylinder C4(p1,p4,rad*0.1,ccol);output << C4.GetCode();
 POVCylinder C5(p1,p5,rad*0.1,ccol);output << C5.GetCode();
 POVCylinder C6(p2,p4,rad*0.1,ccol);output << C6.GetCode();
 POVCylinder C7(p2,p6,rad*0.1,ccol);output << C7.GetCode();
 POVCylinder C8(p3,p5,rad*0.1,ccol);output << C8.GetCode();
 POVCylinder C9(p3,p6,rad*0.1,ccol);output << C9.GetCode();
 POVCylinder C10(p4,p7,rad*0.1,ccol);output << C10.GetCode();
 POVCylinder C11(p5,p7,rad*0.1,ccol);output << C11.GetCode();
 POVCylinder C12(p6,p7,rad*0.1,ccol);output << C12.GetCode();
 //Text
 //FIXME : Corregir settings para valores de texto.
 Vector textpos = b-1*c;
 Vector temppos = 0.25*a+b-1*c; 
 Vector rotate(90-angx*180/M_PI,-(90-angy*180/M_PI),0);
 std::ostringstream text;
 text << "Temperature of system is : " <<sc.Temperature() << " [K]";
 POVText tempe(text.str(),temppos,rotate,0,0,ccol);
 POVText texto("http://www.gnm.cl",textpos,rotate,0,0,ccol);
 output << texto.GetCode();
 output << tempe.GetCode();
 //Set atoms.
 if(level==0)
 {
  for(int i=0;i<sc.Size();i++)
  {
   Vector pos = sc.GetAtom(i).Position();
   rad = ElemRad[sc.GetAtom(i).Species()];
   double mas = ElemMass[sc.GetAtom(i).Species()];
   Vector col(mas,1/rad+1/mas,rad);
   POVSphere S(pos,rad,col);
   output << S.GetCode();
  }
 }
 else
 {
  std::cerr << "NIVEL NO IMPLEMENTADO!!!" << std::endl;
 }
 //close File.
 output.close();
 counter++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new POVRAY(args); }
void destroy(Module * m) { delete m; }


