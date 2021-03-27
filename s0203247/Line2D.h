//
// Created by User on 19/02/2021.
//

#ifndef PROJECT_LINE_H
#define PROJECT_LINE_H
#include "Color.h"
#include "Point2D.h"
#include "easy_image.h"

class Line2D {
    Point2D p1;
    Point2D p2;
    Color color;
public:
    Line2D(Point2D p1, Point2D p2, Color color): p1(p1), p2(p2), color(color){}
    Point2D getP1() const {return p1;}
    Point2D getP2() const {return p2;}
    Color getColor(){return color;}
    img::Color getImageColor() const {return color.imageColor();}

    void reScale(double factor) {
        p1.reScale(factor);
        p2.reScale(factor);
    }

    void move(double dx, double dy) {
        p1.move(dx,dy);
        p2.move(dx,dy);
    }
};


#endif //PROJECT_LINE_H
