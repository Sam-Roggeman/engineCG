//
// Created by User on 19/02/2021.
//

#ifndef PROJECT_POINT_H
#define PROJECT_POINT_H


class Point2D {
    double x;
    double y;

    double z;
public:
    Point2D(double x, double y): x(x), y(y){}
    double getY() const {return y;}
    double getX() const{return x;}

    double getZ() const {
        return z;
    }

    void setZ(double z) {
        Point2D::z = z;
    }

    void reScale(double d) {
        x *=d;
        y *=d;
    }

    void move(double dx, double dy) {
        x += dx;
        y += dy;
    }
};


#endif //PROJECT_POINT_H
