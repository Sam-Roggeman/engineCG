//
// Created by User on 27/03/2021.
//

#include "Line2D.h"
#include "vector3d.h"
#include "matrices.h"
#include <vector>
#include <cmath>
#include "Figuur.h"
using Lines2D = std::list<Line2D>;
using Figures3D = std::list<Figuur>;

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

//void sortp( std::vector<Vector3D>& v){
//    std::vector<Vector3D> sorted = {};
//    Vector3D closest;
//    for (auto it = v.rbegin(); it != v.rend();it++){
//
//    }
//}

void calculateScaleOffset(const Lines2D &lines, const int size, double &d, double& Imagex, double& Imagey, double &dx, double &dy){
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();
    findExtrema(xmin, xmax, ymin, ymax, lines);
    double xrange = xmax-xmin;
    double yrange = ymax-ymin;
    Imagex = size*((xrange)/(fmax(xrange,yrange)));
    Imagey = size*((yrange)/(fmax(xrange,yrange)));
    d = 0.95*(Imagex/xrange);

    double DCx = d*((xmin+xmax)/2);
    double DCy = d*((ymin+ymax)/2);
    dx = Imagex/2 - DCx;
    dy = Imagey/2 - DCy;

}

void generateFractal(Figuur& fig, Figures3D& fractal, const int nr_iterations, const double scale, int count){
    if (count == nr_iterations){
        fractal.emplace_back(fig);
        return;
    }
    Matrix schal_mat = scalingMatrix(1/scale);
    Figuur f_cp = Figuur(fig);
    Matrix trans_mat;
    applyTransformation(f_cp, schal_mat);
    for (int i = 0; i != fig.points.size(); i++){
        Figuur f_cp_cp = Figuur(f_cp);
        trans_mat = translate(fig.points[i]-f_cp_cp.points[i]);
        applyTransformation(f_cp_cp,trans_mat);
        generateFractal(f_cp_cp,fractal,nr_iterations,scale,count+1);
    }
}

void generateFractal(Figuur& fig, Figures3D& fractal, const int nr_iterations, const double scale){
    generateFractal(fig, fractal,nr_iterations,scale, 0);
}


#endif //ENGINE_helpfunctions_H