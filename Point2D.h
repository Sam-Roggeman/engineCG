//
// Created by User on 19/02/2021.
//

#ifndef PROJECT_POINT_H
#define PROJECT_POINT_H


class Point2D {
    double x;
    double y;
public:
    Point2D(double x, double y): x(x), y(y){}
    double getY() const {return y;}
    double getX() const{return x;}
};


#endif //PROJECT_POINT_H
