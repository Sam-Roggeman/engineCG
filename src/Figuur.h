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
    Figuur()= default;;
    explicit Figuur(Color color) : ambientReflection(color) {};
    std::vector<Vector3D> points;
    std::vector<Vlak> vlakken;
//    Color color = Color(0, 0, 0);
    Color ambientReflection = Color(0,0,0);
    Color diffuseReflection = Color(0,0,0);
    Color specularReflection = Color(0,0,0);
    double reflectionCoefficient = 1;

    bool hasAmb(){
        return !ambientReflection.isZero();
    }

    bool hasDiff(){
        return !diffuseReflection.isZero();
    }

    bool hasSpec(){
        return !specularReflection.isZero();
    }


    void clear(){
        points.clear();
        vlakken.clear();
    }


    void addpoint(double x, double y, double z) {
        auto p = Vector3D(x,y,z,false);
        points.emplace_back(p);
    }

    Figuur(const Color &ambientReflection, const Color &diffuseReflection, const Color &specularReflection,
           double &reflectionCoefficient) : ambientReflection(ambientReflection), diffuseReflection(diffuseReflection),
                                           specularReflection(specularReflection),
                                           reflectionCoefficient(reflectionCoefficient) {}

    Figuur(const Figuur& f){
        this->clear();
        this->points.insert(this->points.end(),f.points.begin(),f.points.end());
        this->vlakken.insert(this->vlakken.end(),f.vlakken.begin(),f.vlakken.end());

        ambientReflection = f.ambientReflection;
        specularReflection = f.specularReflection;
        diffuseReflection = f.diffuseReflection;
        reflectionCoefficient = f.reflectionCoefficient;

    }

    void addfigure(const Figuur& figuur2) {
        int current_point_ind = points.size();
        std::vector<int> vlak_ind;
        this->points.insert(this->points.end(),figuur2.points.begin(),figuur2.points.end());
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

    void addVlak(const std::vector<Vector3D> &vector) {
        std::vector<int> ind;
        for (auto& hoek:vector){
            ind.emplace_back(points.size());
            points.emplace_back(hoek);
        }
        vlakken.emplace_back(Vlak(ind));
    }

    void triangulate()  {
        auto temp = vlakken;
        vlakken.clear();
        Vlak n_vlak;
        int size;
        for (auto &vlak:temp) {
            size = vlak.point_indexes.size();
            for (int ind = 1; ind < size-1; ind++){
                n_vlak = Vlak({vlak.point_indexes[0], vlak.point_indexes[ind], vlak.point_indexes[ind+1]});
                vlakken.emplace_back(n_vlak);
            }
        }

    }

};


#endif //ENGINE_FIGUUR_H
