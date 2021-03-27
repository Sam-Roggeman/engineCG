//
// Created by User on 27/03/2021.
//

#include "Line2D.h"
#include <vector>
#include "Line2D.h"
#include <cmath>
using Lines2D = std::list<Line2D>;

#ifndef ENGINE_helpfunctions_H
#define ENGINE_helpfunctions_H

inline int roundToInt(double d)
{
    return static_cast<int>(std::round(d));
}


void findExtrema(double &xmin, double &xmax, double &ymin, double &ymax, const Lines2D& lines){
    if ( lines.empty()){std::cerr<<"No lines" << std::endl;}
    for (const auto &line:lines){
        Point2D p = line.getP1();
        for (int i = 0; i<2;i++){
            if (p.getX() < xmin) xmin = p.getX();
            else if (p.getX() > xmax) xmax = p.getX();
            if (p.getY() < ymin) ymin = p.getY();
            else if (p.getY() > ymax) ymax = p.getY();
            p = line.getP2();
        }
    }
}

#endif //ENGINE_helpfunctions_H