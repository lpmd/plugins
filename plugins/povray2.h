//
//
//

#ifndef __POVRAY2_H__
#define __POVRAY2_H__

#include <lpmd/visualizer.h>
#include <lpmd/plugin.h>
#include <map>

using namespace lpmd;

class POVRAY2: public lpmd::Visualizer, public lpmd::Module
{
 public:
    POVRAY2(std::string args);
    ~POVRAY2();

    // From Module
    void SetParameter(std::string name);
    void Show() const;
    void ShowHelp() const;
    std::string Keywords() const;

    // From Visualizer
    void Apply(const lpmd::MD & md);

 private:
    //Informacion Directorios
    std::string header;
    std::string direct;
    //Informacion cuadros de Texto.
    std::string mtexts[20];
    std::string pos[20];
    std::string colors[20];
    double scale[20];
    std::string extras[20];
    int ntext;
    //Informacion del logo
    std::string file_logo;
    std::string pos_logo;
    double scale_logo;
    //Informacion Celda de Simulacion.
    std::string angle_cell;
    bool show_cell;
    std::string type_cell;
    std::string color_cell;
    double scale_cell;
    //General
    bool campos;
    Vector original_camera;
    Vector original_lok;
    Vector original_typetext;
    std::string background;
    //Contador de Ficheros POV y Nivel Especial.
    int counter;
    int level;
};

#endif

