//
// Created by User on 27/03/2021.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H
#include <vector>
#include <limits>
#include <iostream>

class ZBuffer: public std::vector<std::vector<double>> {
public:

    ZBuffer(const int width, const int height){
        (*this).resize(height);
        for (std::vector<double> &hor:*this) {
            hor.resize(width);
            for (double &pixel_z:hor) {
                pixel_z = std::numeric_limits<double>::max();
            }
        }
    }

    bool changeIfCloser(unsigned int row, unsigned int column, double z_inverse){

        double *current = &((*this).at(row).at(column));
        if (*current>z_inverse) {
            *current = z_inverse;
            return true;
        }
        return false;
    };
};


#endif //ENGINE_ZBUFFER_H
