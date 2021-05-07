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
    Color color = Color(0, 0, 0);

    void addpoint(double x, double y, double z) {
        auto p = Vector3D(x,y,z,false);
        points.emplace_back(p);
    }

    Figuur(const Figuur& f){
        for (const auto &point: f.points){
            this->points.emplace_back(Vector3D(point));
        }
        for (const auto& vlak:f.vlakken ){
            vlakken.emplace_back(Vlak(vlak));
        }
        color = f.color;

    }

    void addfigure(const Figuur& figuur2) {
        int current_point_ind = points.size();
        std::vector<int> vlak_ind;
        for (const auto& point : figuur2.points){
            points.emplace_back(point);
        }
        for (const auto& vlak:figuur2.vlakken){
            vlak_ind = {};
            for (const auto& point_ind: vlak.point_indexes){
                vlak_ind.emplace_back(current_point_ind+point_ind);
            }

            vlakken.emplace_back(vlak_ind);
        }
    }

    std::vector<Vector3D> vlakInPoints(int index)const{
        std::vector<Vector3D> v;
        for (int p_ind: vlakken[index].point_indexes){
            v.emplace_back(points[p_ind]);
        }
        return v;
    }

    void addVlak(std::vector<Vector3D> &vector) {
        std::vector<int> ind;
        for (auto& hoek:vector){
            ind.emplace_back(points.size());
            points.emplace_back(hoek);
        }
        vlakken.emplace_back(Vlak(ind));

    }
};


#endif //ENGINE_FIGUUR_H
