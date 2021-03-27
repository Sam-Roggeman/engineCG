//
// Created by User on 27/03/2021.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H
#include <vector>
#include <limits>

class ZBuffer: public std::vector<std::vector<double>> {

    const int height;
    const int width;
public:

    ZBuffer(const int width, const int height) : width(width), height(height){
        this->reserve(height);
        for (std::vector<double> &hor:*this){
            hor.reserve(width);
            for (double &pixel_z:hor){
                pixel_z = std::numeric_limits<double>::max();
            }
    }
};


#endif //ENGINE_ZBUFFER_H
