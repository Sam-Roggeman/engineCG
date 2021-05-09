//
// Created by User on 19/02/2021.
//

#ifndef PROJECT_COLOR_H
#define PROJECT_COLOR_H
#include "easy_image.h"

class Color {
    double red;
    double green;
    double blue;
public:
    Color(double red,  double blue,double green): red(red),  green(green),blue(blue){};
    explicit Color(std::vector<double> rgb){red = rgb[0]; green=rgb[1]; blue = rgb[2]; }
    img::Color imageColor() const{return img::Color(red*255, green*255, blue*255);}
    Color operator*(const Color &c)const{
        return Color(red*c.red, blue*c.blue, green*c.green);
    }
    Color operator*(const double d)const{
        return Color(red*d, blue*d, green*d);
    }
    Color& operator+=(Color c){
        red = std::min(red+c.red,1.0);
        blue = std::min(blue+c.blue,1.0);
        green = std::min(green+c.green,1.0);
        return *this;
    }
};


#endif //PROJECT_COLOR_H
