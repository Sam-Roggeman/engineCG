//
// Created by User on 8/03/2021.
//

#ifndef ENGINE_FIGUUR_H
#define ENGINE_FIGUUR_H
#include <vector>
#include "vector3d.h"
#include "Color.h"
#include "Vlak.h"

class Figuur {
public:
    explicit Figuur(Color color) : color(color) {};
    std::vector<Vector3D> points;
    std::vector<Vlak> vlakken;
    Color color;

    void addpoint(double x, double y, double z) {
        Vector3D* p = new Vector3D(x,y,z,false);
        points.emplace_back(*p);
    }
};


#endif //ENGINE_FIGUUR_H
