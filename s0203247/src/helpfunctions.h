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

std::vector<Vlak> triangulate(const Vlak& v)  {
    std::vector<Vlak> vlakken = {};
    Vlak n_vlak;
    int size = v.point_indexes.size();

    for (int ind = 1; ind < size-1; ind++){
        n_vlak = Vlak({v.point_indexes[0], v.point_indexes[ind], v.point_indexes[ind+1]});
        vlakken.emplace_back(n_vlak);
    }
    return vlakken;
}

//void findExtrema(double &xmin, double &xmax, double &ymin, double &ymax, const std::vector<Figuur> figuren){
//    if ( figuren.empty()){std::cerr<<"No lines" << std::endl;}
//    for (const auto &figuur:figuren){
//        Point2D p = figuur.;
//        for (int i = 0; i<2;i++){
//            if (p.getX() < xmin) xmin = p.getX();
//            else if (p.getX() > xmax) xmax = p.getX();
//            if (p.getY() < ymin) ymin = p.getY();
//            else if (p.getY() > ymax) ymax = p.getY();
//            p = line.getP2();
//        }
//    }
//}

void findYExtrema(double &ymin, double &ymax, const Lines2D& lines){
    Point2D p = lines.front().getP1();
    if ( lines.empty()){std::cerr<<"No lines" << std::endl;}
    for (const auto &line:lines){
        p = line.getP1();
        for (int i = 0; i<2;i++){
            if (p.getY() < ymin) ymin = p.getY();
            else if (p.getY() > ymax) ymax = p.getY();
            p = line.getP2();
        }
    }
}

#endif //ENGINE_helpfunctions_H