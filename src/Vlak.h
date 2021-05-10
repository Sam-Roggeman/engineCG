//
// Created by User on 8/03/2021.
//

#ifndef ENGINE_VLAK_H
#define ENGINE_VLAK_H
#include <utility>
#include <vector>

class Vlak {
public:
    std::vector<int> point_indexes;

    Vlak(){point_indexes={};};
    explicit Vlak(std::vector<int> point_indexes) : point_indexes(std::move(point_indexes)){};
    Vlak(const Vlak& vlak){
        for (const auto& point_ind : vlak.point_indexes) {
            point_indexes.emplace_back(point_ind);
        }
    }
};


#endif //ENGINE_VLAK_H
