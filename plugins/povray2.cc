//
//
//


#include <sstream>

#include <lpmd/util.h>
#include <lpmd/inputfile.h>

#include "povray2.h"

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
  POVText(std::string t,Vector p,Vector r,double x,double y, Vector d,double s){text=t;pos=p;rot=r;a=x,b=y;col=d;scale=s;}
  ~POVText(){};
  std::string GetCode()
  {
   std::ostringstream ostr;
   ostr << "text { \n";
   ostr << "\t ttf\"timrom.ttf\" " << text << " "<<a<<", "<<b<<"\n";
   ostr << "\t pigment {color "<<col.GetPOV()<<"}\n";
   ostr << "\t scale     "<<scale<<"\n";
   ostr << "\t rotate    "<<rot.GetPOV()<<"\n";
   ostr << "\t translate "<<pos.GetPOV()<<"\n";
   ostr << "}" << '\n';
   return ostr.str();
  }
 private:
  Vector pos,col,rot;
  double a,b,scale;
  std::string text;
};

class POVLogo
{
 public:
  POVLogo(std::string f,Vector p,double s) {file = f; pos = p; scale = s;}
  ~POVLogo(){};
  std::string GetCode()
  {
   std::ostringstream ostr;
   ostr << "#declare logo=box { <0,0,0> <1,1,0> \n";
   //ostr << "texture { Chrome_Texture }\n";
   ostr << "texture { \n";
   ostr << "    pigment { \n";
   ostr << "        image_map { \n";
   ostr << "             gif " << file << "\n";
   ostr << "             map_type 0 \n";
   ostr << "             interpolate 2 \n";
   ostr << "             once \n";
   ostr << "                  } \n";
   ostr << "            } \n";
   //ostr << "        finish {diffuse 0.9}\n";
   //ostr << "        translate 0.1*x+0.1*y \n";
   ostr << "        } \n";
   ostr << "scale " << scale << " \n";
   ostr << "texture {Glass} \n";
   ostr << "}\n";
   ostr << "object {logo translate " << pos.GetPOV() << "  } \n";
   return ostr.str();
  }
 private:
  std::string file;
  Vector pos;
  double scale;
};

POVRAY2::POVRAY2(std::string args): Module("povray")
{
 direct = "povray";
 header = "movie-";
 counter = 0;
 level=0;
 ntext=0;
 background = "<0,0,0>" ; 
 angle_cell = "<-20,-20,0>" ;
 pos_logo = "<ur>" ;
 scale_logo = 1.0 ;
 file_logo = "";
 campos=true;
 ProcessArguments(args);
}

POVRAY2::~POVRAY2()
{
}

void POVRAY2::SetParameter(std::string name)
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
 if (name == "text")
 {
  std::string texto = GetNextWord('"');
  std::string posic = GetNextWord('<');
  std::string color = GetNextWord('<');
  std::string scal = GetNextWord('[');
  std::string extra = GetNextWord('(');
  mtexts[ntext]=texto;
  pos[ntext]=posic;
  colors[ntext]=color;
  extras[ntext]=extra;
  //Borra los bordes, solo scale y extra.
  scal.erase(0,1);
  scal.erase(scal.size()-1,1);
  scale[ntext] = atof(scal.c_str());
  extras[ntext].erase(0,1);
  extras[ntext].erase(extras[ntext].size()-1,1);
  ntext++;
 }
 if (name == "logo")
 {
  AssignParameter("logo-file",GetNextWord('"'));
  AssignParameter("logo-scale",GetNextWord());
  AssignParameter("logo-pos",GetNextWord('<'));
  file_logo = GetString("logo-file");
  scale_logo = GetDouble("logo-scale");
  pos_logo = GetString("logo-pos");
 }
 if (name == "rotate")
 {
  AssignParameter("rotate",GetNextWord('<'));
  angle_cell = GetString("rotate");
 }
 if (name == "box")
 {
  AssignParameter("cell-show",GetNextWord());
  if(GetBool("cell-show")!=false)
  {
  AssignParameter("cell-type",GetNextWord());
  AssignParameter("cell-color",GetNextWord('<'));
  AssignParameter("cell-scale",GetNextWord());
  show_cell = GetBool("cell-show");
  type_cell = GetString("cell-type");
  color_cell = GetString("cell-color");
  scale_cell = GetDouble("cell-scale");
  }
 }
 if (name == "background")
 {
  AssignParameter("background",GetNextWord('<'));
  background = GetString("background");
 }
 if (name == "camera")
 {
  AssignParameter("camera",GetNextWord());
  campos = GetBool("camera");
 }
}

void POVRAY2::Show() const
{
 Module::Show();
 std::cout << "   start    = " << start_step << '\n';
 std::cout << "   end      = " << end_step << '\n';
 std::cout << "   step     = " << interval << '\n';
 std::cout << "   header   = " << header << '\n';
 std::cout << "   direct   = " << direct << '\n';
}

void POVRAY2::ShowHelp() const
{
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << " Module Name        = POVRAY2 " << '\n';
 std::cout << " Module Version     = 1.0 " << '\n';
 std::cout << " Support API lpmd   = 1.0 " << '\n';
 std::cout << " Problems Report to = gnm@gnm.cl " << '\n';
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
 std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

std::string POVRAY2::Keywords() const { return "start end step header direct texti logo rotate box background camera"; }

void POVRAY2::Apply(const MD & md)
{
 //Carga celda de simulacion y genera directorio/archivo para guardar ficheros .pov
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
 //Set de elemtons basicos de cada fichero .pov y variables globales del sistema.
 output << "#include \"colors.inc\"" << '\n';
 output << "#include \"stones.inc\"" << '\n';
 output << "#include \"textures.inc\"" << '\n';
 output << "#include \"shapes.inc\"" << '\n';
 output << "#include \"glass.inc\"" << '\n';
 output << "#include \"metals.inc\"" << '\n';
 output << "#include \"woods.inc\"" << '\n';
 Vector red(1,0,0);
 Vector green(0,1,0);
 Vector blue(0,0,1);
 //Set Luces y Camara
 Vector ratec(0,0,-1*(b.Mod()/c.Mod()+a.Mod()/c.Mod()));
 Vector rateb(0,-1*(c.Mod()/b.Mod()),0);
 Vector ratea(-1*(c.Mod()/a.Mod()),0,0);
 Vector camerapos(0,0,-1*(b.Mod()/c.Mod()+a.Mod()/c.Mod())*c.Mod());
 Vector translate = (-0.5*a-0.5*b+1*c);
 Vector cameralok(0,0,c.Mod());
 Vector corner1(0,0,0);
 Vector corner2 = a+b+c;
 Vector light1=Vector(camerapos.GetX(),camerapos.GetY(),camerapos.GetZ());
 Vector typetext(a.Mod()*ratea.Mod(),b.Mod()*rateb.Mod(),-(camerapos+cameralok).Mod());
 if(counter == 0 && campos == false)
 {
  original_camera = camerapos;
  original_lok = cameralok;
  original_typetext = typetext;
 }
 if(counter > 0 && campos == false) {camerapos = original_camera; cameralok=original_lok;typetext = original_typetext;}
 light1=Vector(camerapos.GetX(),camerapos.GetY(),camerapos.GetZ());
 output << "camera {\n";
 output << "\t location " << camerapos.GetPOV() << '\n';
 output << "\t look_at  " << cameralok.GetPOV() << '\n';
 output << "}\n\n";
 output << "light_source { " << light1.GetPOV() << "color White}"<< '\n';
 output << "background   { color " << background << " } " << '\n'; 
 //Set de Textos ingresados en el modulo
 Vector ur = typetext*Vector( 0.4, 0.4, 1);
 Vector uc = typetext*Vector(-0.1, 0.4, 1);
 Vector ul = typetext*Vector(-0.6, 0.4, 1);
 Vector cl = typetext*Vector(-0.6, 0.0, 1);
 Vector cr = typetext*Vector( 0.4, 0.0, 1);
 Vector dr = typetext*Vector( 0.4,-0.4, 1);
 Vector dc = typetext*Vector(-0.1,-0.4, 1);
 Vector dl = typetext*Vector(-0.6,-0.4, 1);
 for(int i=0;i<ntext;i++)
 {
  Vector tp,tc;
  //POSICION
  if(pos[i]=="<ur>") tp=ur;
  else if(pos[i]=="<uc>") tp=uc;
  else if(pos[i]=="<ul>") tp=ul;
  else if(pos[i]=="<cl>") tp=cl;
  else if(pos[i]=="<cr>") tp=cr;
  else if(pos[i]=="<dr>") tp=dr;
  else if(pos[i]=="<dl>") tp=dl;
  else if(pos[i]=="<dc>") tp=dc;
  else {tp = Vector(pos[i]);}
  //COLOR
  if(colors[i]=="<red>") tc=Vector(1,0,0);
  else if(colors[i]=="<blue>") tc=Vector(0,0,1);
  else if(colors[i]=="<green>") tc=Vector(0,1,0);
  else {tc = Vector(colors[i]);}
  //SCALE
  double sca = scale[i];
  //EXTRAS
  std::string phrase=mtexts[i];
  std::vector<std::string> EV,TV;
  if(extras[i].size()>0) EV = SplitTextLine(extras[i],',');
  if(mtexts[i].size()>0) TV = SplitTextLine(mtexts[i],'%');
  ostr.str("");

  if (EV.size()==TV.size()-1 && EV.size()>0)
  {
   double dval=0.0e0;
   Vector vval(0,0,0);
   for(unsigned int k=0;k<EV.size();k++)
   {
    if (EV[k] == "Temp" || EV[k] == "temp") dval=sc.Temperature();
    else if(EV[k] == "Vol" || EV[k] == "vol") dval=sc.Volume();
    else if(EV[k] == "Step" || EV[k] == "step") dval = md.CurrentStep();
    else {dval=0.0e0;}

    if(k==0) ostr << TV[k] <<" "<< dval << " "<< TV[k+1];
    else {ostr << " " << dval << " " << TV[k+1];}
    phrase=ostr.str();
   }
  }
  if(EV.size()>0 && TV.size()<EV.size()+1) {std::cerr << "WARNING : POVRAY MODULE HAVE MORE VARIBLAES THAT % SYMBOLS\n";}
  //CONSTRUCCION
  POVText t(phrase,tp,Vector(0,0,0),0.1,0.0,tc,sca);
  output << t.GetCode();
 }
 //LOGO (ES UN OBJECT POR CONSTRUCCION).
 if(file_logo!="")
 {
  Vector flp;
  if(pos_logo=="<ur>") flp=ur;
  else if(pos_logo=="<uc>") flp=uc;
  else if(pos_logo=="<ul>") flp=ul;
  else if(pos_logo=="<cl>") flp=cl;
  else if(pos_logo=="<cr>") flp=cr;
  else if(pos_logo=="<dr>") flp=dr;
  else if(pos_logo=="<dl>") flp=dl;
  else if(pos_logo=="<dc>") flp=dc;
  else {flp = Vector(pos_logo);}
  POVLogo logo(file_logo,flp,scale_logo);
  output << logo.GetCode();
 }
 //CELDA DE SIMULACION COMO UNA UNION
 output << "#declare cell = union {\n";
 //Set Cilindros
 Vector p0(0,0,0);
 Vector p1=a,p2=b,p3=c;
 Vector p4=p1+p2,p5=p1+p3,p6=p3+p2,p7=p1+p2+p3;
 double rad = ElemRad[sc.GetAtom(0).Species()];
 POVCylinder C1(p0,p1,rad*0.1,red);output << C1.GetCode();
 POVCylinder C2(p0,p2,rad*0.1,red);output << C2.GetCode();
 POVCylinder C3(p0,p3,rad*0.1,red);output << C3.GetCode();
 POVCylinder C4(p1,p4,rad*0.1,red);output << C4.GetCode();
 POVCylinder C5(p1,p5,rad*0.1,red);output << C5.GetCode();
 POVCylinder C6(p2,p4,rad*0.1,red);output << C6.GetCode();
 POVCylinder C7(p2,p6,rad*0.1,red);output << C7.GetCode();
 POVCylinder C8(p3,p5,rad*0.1,red);output << C8.GetCode();
 POVCylinder C9(p3,p6,rad*0.1,red);output << C9.GetCode();
 POVCylinder C10(p4,p7,rad*0.1,red);output << C10.GetCode();
 POVCylinder C11(p5,p7,rad*0.1,red);output << C11.GetCode();
 POVCylinder C12(p6,p7,rad*0.1,red);output << C12.GetCode();
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
 output << "}";
 output << "object { cell ";
 output << "rotate "<<angle_cell<<'\n';
 output << "translate "<<translate.GetPOV() << "}";
 //close File.
 output.close();
 counter++;
}

// Esto se incluye para que el modulo pueda ser cargado dinamicamente
Module * create(std::string args) { return new POVRAY2(args); }
void destroy(Module * m) { delete m; }


