//
// Created by User on 27/03/2021.
//

#ifndef ENGINE_PREDMADEFIGURES_H
#define ENGINE_PREDMADEFIGURES_H
#include "matrices.h"
#include <unordered_map>
#include "helpfunctions.h"
Figuur createCube(Color c){
    Figuur f = Figuur(c);
    f.points.emplace_back(Vector3D::point(1,-1,-1));
    f.points.emplace_back(Vector3D::point(-1,1,-1));
    f.points.emplace_back(Vector3D::point(1,1,1));
    f.points.emplace_back(Vector3D::point(-1,-1,1));
    f.points.emplace_back(Vector3D::point(1,1,-1));
    f.points.emplace_back(Vector3D::point(-1,-1,-1));
    f.points.emplace_back(Vector3D::point(1,-1,1));
    f.points.emplace_back(Vector3D::point(-1,1,1));
    f.vlakken.emplace_back(Vlak({0,4,2,6}));
    f.vlakken.emplace_back(Vlak({4,1,7,2}));
    f.vlakken.emplace_back(Vlak({1,5,3,7}));
    f.vlakken.emplace_back(Vlak({5,0,6,3}));
    f.vlakken.emplace_back(Vlak({6,2,7,3}));
    f.vlakken.emplace_back(Vlak({0,5,1,4}));
    return f;
}

Figuur createTetrahedron(Color color) {
    Figuur f = Figuur(color);
    f.points.emplace_back(Vector3D::point(1,-1,-1));
    f.points.emplace_back(Vector3D::point(-1,1,-1));
    f.points.emplace_back(Vector3D::point(1,1,1));
    f.points.emplace_back(Vector3D::point(-1,-1,1));
    f.vlakken.emplace_back(Vlak({0,1,2}));
    f.vlakken.emplace_back(Vlak({1,3,2}));
    f.vlakken.emplace_back(Vlak({0,3,1}));
    f.vlakken.emplace_back(Vlak({0,2,3}));
    return f;
}

Figuur createOctahedron(Color color) {
    Figuur f = Figuur(color);
    f.points.emplace_back(Vector3D::point(1,0,0));
    f.points.emplace_back(Vector3D::point(0,1,0));
    f.points.emplace_back(Vector3D::point(-1,0,0));
    f.points.emplace_back(Vector3D::point(0,-1,0));
    f.points.emplace_back(Vector3D::point(0,0,-1));
    f.points.emplace_back(Vector3D::point(0,0,1));
    f.vlakken.emplace_back(Vlak({0,1,5}));
    f.vlakken.emplace_back(Vlak({1,2,5}));
    f.vlakken.emplace_back(Vlak({2,3,5}));
    f.vlakken.emplace_back(Vlak({3,0,5}));
    f.vlakken.emplace_back(Vlak({1,0,4}));
    f.vlakken.emplace_back(Vlak({2,1,4}));
    f.vlakken.emplace_back(Vlak({3,2,4}));
    f.vlakken.emplace_back(Vlak({0,3,4}));
    return f;
}

Figuur createIcosahedron(Color color) {
    Figuur f = Figuur(color);
    f.points.emplace_back(Vector3D::point(0,0,sqrt(5)/2));
    for (int i = 2; i< 7;i++){
        double d = (i-2)*2*M_PI/5;
        f.points.emplace_back(Vector3D::point(cos(d), sin(d),0.5));
    }
    for (int i = 7; i< 12;i++){
        double d = (i-7)*2*M_PI/5+M_PI/5;
        f.points.emplace_back(Vector3D::point(cos(d), sin(d),-0.5));
    }
    f.points.emplace_back(Vector3D::point(0,0,-sqrt(5)/2));

    f.vlakken.emplace_back(Vlak({0,1,2}));
    f.vlakken.emplace_back(Vlak({0,2,3}));
    f.vlakken.emplace_back(Vlak({0,3,4}));
    f.vlakken.emplace_back(Vlak({0,4,5}));
    f.vlakken.emplace_back(Vlak({0,5,1}));

    f.vlakken.emplace_back(Vlak({1,6,2}));
    f.vlakken.emplace_back(Vlak({2,6,7}));
    f.vlakken.emplace_back(Vlak({2,7,3}));
    f.vlakken.emplace_back(Vlak({3,7,8}));
    f.vlakken.emplace_back(Vlak({3,8,4}));

    f.vlakken.emplace_back(Vlak({4,8,9}));
    f.vlakken.emplace_back(Vlak({4,9,5}));
    f.vlakken.emplace_back(Vlak({5,9,10}));
    f.vlakken.emplace_back(Vlak({5,10,1}));
    f.vlakken.emplace_back(Vlak({1,10,6}));

    f.vlakken.emplace_back(Vlak({11,7,6}));
    f.vlakken.emplace_back(Vlak({11,8,7}));
    f.vlakken.emplace_back(Vlak({11,9,8}));
    f.vlakken.emplace_back(Vlak({11,10,9}));
    f.vlakken.emplace_back(Vlak({11,6,10}));

    return f;
}

Figuur createBuckyball(Color color) {
    Figuur f = createIcosahedron(color);
    Vector3D result;
    Vector3D result2;
    Vector3D start;
    Vector3D end;
    const Vector3D* start_ptr;
    std::unordered_map<const Vector3D*, std::vector<Vector3D>> temppoints;
    std::vector<int> temp;
    std::vector<int> temp2;
    Figuur bucky = Figuur(color);
    for (auto & vlak:f.vlakken){
        temp = {};
        for (auto p_ind = 0; p_ind < vlak.point_indexes.size(); p_ind++){
            start = f.points[vlak.point_indexes[p_ind]];
            end = f.points[vlak.point_indexes[(p_ind+1)%vlak.point_indexes.size()]];
            temp.emplace_back(bucky.points.size());
            temp.emplace_back(bucky.points.size()+1);
            bucky.points.emplace_back(Vector3D(start.x+(end.x-start.x)/3,start.y+(end.y-start.y)/3,start.z+(end.z-start.z)/3,false));
            bucky.points.emplace_back(Vector3D(start.x+(end.x-start.x)*2/3,start.y+(end.y-start.y)*2/3,start.z+(end.z-start.z)*2/3,false));
        }
        bucky.vlakken.emplace_back(Vlak(temp));
        temp = {};
        for (auto p_ind = 0; p_ind < vlak.point_indexes.size(); p_ind++){
            temp = {};
            start_ptr = &f.points[vlak.point_indexes[p_ind]];
            end = f.points[vlak.point_indexes[(p_ind+1)%vlak.point_indexes.size()]];
            temppoints[start_ptr].emplace_back(Vector3D(start_ptr->x+(end.x-start_ptr->x)/3,start_ptr->y+(end.y-start_ptr->y)/3,start_ptr->z+(end.z-start_ptr->z)/3,false)) ;
        }
    }
    for (const auto& point : f.points) {
        for (auto &buckypoint : temppoints.at(&point)) {
            bucky.points.emplace_back(buckypoint);
        }
    }
    bucky.vlakken.emplace_back(Vlak({120,121,122,123,124}));
    bucky.vlakken.emplace_back(Vlak({125,126,128,129,127}));
    bucky.vlakken.emplace_back(Vlak({130,131,134,133,132}));
    bucky.vlakken.emplace_back(Vlak({135,136,139,138,137}));
    bucky.vlakken.emplace_back(Vlak({140,141,144,143,142}));
    bucky.vlakken.emplace_back(Vlak({145,146,149,148,147}));
    bucky.vlakken.emplace_back(Vlak({150,151,153,154,152}));
    bucky.vlakken.emplace_back(Vlak({155,156,157,159,158}));
    bucky.vlakken.emplace_back(Vlak({160,161,162,164,163}));
    bucky.vlakken.emplace_back(Vlak({165,166,167,169,168}));
    bucky.vlakken.emplace_back(Vlak({170,171,172,174,173}));
    bucky.vlakken.emplace_back(Vlak({175,176,177,178,179}));

    return bucky;
}

Figuur createDodecahedron(Color color) {
    Figuur f = Figuur(color);
    Figuur ico = createIcosahedron(color);
    for (const Vlak& ico_vlak:ico.vlakken){
        double x = 0;
        double y = 0;
        double z = 0;
        for (int ico_ind:ico_vlak.point_indexes) {
            x += ico.points[ico_ind].x;
            y += ico.points[ico_ind].y;
            z += ico.points[ico_ind].z;
        }
        f.points.emplace_back(Vector3D::point(x/3,y/3,z/3));
    }
    f.vlakken.emplace_back(Vlak({0,1,2,3,4}));
    f.vlakken.emplace_back(Vlak({0,5,6,7,1}));
    f.vlakken.emplace_back(Vlak({1,7,8,9,2}));
    f.vlakken.emplace_back(Vlak({2,9,10,11,3}));
    f.vlakken.emplace_back(Vlak({3,11,12,13,4}));

    f.vlakken.emplace_back(Vlak({4,13,14,5,0}));
    f.vlakken.emplace_back(Vlak({19,18,17,16,15}));
    f.vlakken.emplace_back(Vlak({19,14,13,12,18}));
    f.vlakken.emplace_back(Vlak({18,12,11,10,17}));
    f.vlakken.emplace_back(Vlak({17,10,9,8,16}));

    f.vlakken.emplace_back(Vlak({16,8,7,6,15}));
    f.vlakken.emplace_back(Vlak({15,6,5,14,19}));

    return f;
}

Figuur createSphere(int iterations, Color color) {
    Figuur f = createIcosahedron(color);
    int nr = 0;
    while (nr < iterations) {
        std::vector<Vlak> new_vlakken = {};
        nr++;
        for (const Vlak &vlak:f.vlakken) {
            int p_ind1 = vlak.point_indexes[0];
            int p_ind2 = vlak.point_indexes[1];
            int p_ind3 = vlak.point_indexes[2];
            Vector3D p1 = f.points[p_ind1];
            Vector3D p2 = f.points[p_ind2];
            Vector3D p3 = f.points[p_ind3];
            Vector3D p4 = Vector3D::point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2);
            Vector3D p5 = Vector3D::point((p1.x + p3.x) / 2, (p1.y + p3.y) / 2, (p1.z + p3.z) / 2);
            Vector3D p6 = Vector3D::point((p2.x + p3.x) / 2, (p2.y + p3.y) / 2, (p2.z + p3.z) / 2);
            int p_ind4 = f.points.size();
            int p_ind5 = f.points.size()+1;
            int p_ind6 = f.points.size()+2;
            new_vlakken.emplace_back(Vlak({p_ind1,p_ind4,p_ind5}));
            new_vlakken.emplace_back(Vlak({p_ind2,p_ind6,p_ind4}));
            new_vlakken.emplace_back(Vlak({p_ind3,p_ind5,p_ind6}));
            new_vlakken.emplace_back(Vlak({p_ind4,p_ind6,p_ind5}));
            f.points.emplace_back(p4);
            f.points.emplace_back(p5);
            f.points.emplace_back(p6);
        }
        f.vlakken = new_vlakken;
    }
    for (Vector3D& point:f.points){
        point.normalise();
    }
    return f;
}

Figuur createCone(const double height,const int nr_vlakken, Color c){
    Figuur f = Figuur(c);
    for (int i = 0; i < nr_vlakken; i++){
        f.points.emplace_back(Vector3D::point(cos(2*i*M_PI/nr_vlakken),sin(2*i*M_PI/nr_vlakken),0));
    }
    f.points.emplace_back(Vector3D::point(0,0,height));
    std::vector<int> grd_ind = {};
    for (int i = 0; i < nr_vlakken; i++) {
        f.vlakken.emplace_back(Vlak({i, (i+1)%(nr_vlakken),nr_vlakken}));
        grd_ind.emplace_back(i);
    }
    f.vlakken.emplace_back(grd_ind);

    return f;
}

Figuur createCylinder(const double height,const int nr_vlakken, Color c){
    Figuur f = Figuur(c);
    for (int i = 0; i < nr_vlakken; i++){
        f.points.emplace_back(Vector3D::point(cos(2*i*M_PI/nr_vlakken),sin(2*i*M_PI/nr_vlakken),0));
    }
    for (int i = 0; i < nr_vlakken; i++){
        f.points.emplace_back(Vector3D::point(cos(2*i*M_PI/nr_vlakken),sin(2*i*M_PI/nr_vlakken),height));
    }
    std::vector<int> grd_ind = {};
    std::vector<int> bov_ind = {};
    for (int i = 0; i < nr_vlakken; i++) {
        f.vlakken.emplace_back(Vlak({i, (i+1)%(nr_vlakken), nr_vlakken+(i+1)%(nr_vlakken),i+nr_vlakken}));
        grd_ind.emplace_back(i);
        bov_ind.emplace_back(i+nr_vlakken);
    }
    f.vlakken.emplace_back(grd_ind);
    f.vlakken.emplace_back(bov_ind);
    return f;
}

Figuur createTorus(double r, double R, int n, int m, Color color) {
    Figuur f = Figuur(color);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) {
            double u = (2*i*M_PI)/n;
            double v = (2*j*M_PI)/m;
            f.points.emplace_back(Vector3D::point((R+r*cos(v))*cos(u),(R+r*cos(v))*sin(u),r*sin(v)));
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            f.vlakken.emplace_back(Vlak({i*n+j, n*((i+1)%n) + j, n*((i+1)%n) + (j+1)%m ,n*i+(j+1)%m}));
        }
    }
    return f;
}

Figuur createVrachtwagen(Color color){
    Figuur f = Figuur(color);
    f.points.emplace_back(Vector3D::point(1,-2,-1));
    f.points.emplace_back(Vector3D::point(-1,2,-1));
    f.points.emplace_back(Vector3D::point(1,2,1));
    f.points.emplace_back(Vector3D::point(-1,-2,1));
    f.points.emplace_back(Vector3D::point(1,2,-1));
    f.points.emplace_back(Vector3D::point(-1,-2,-1));
    f.points.emplace_back(Vector3D::point(1,-2,1));
    f.points.emplace_back(Vector3D::point(-1,2,1));
    f.vlakken.emplace_back(Vlak({0,4,2,6}));
    f.vlakken.emplace_back(Vlak({4,1,7,2}));
    f.vlakken.emplace_back(Vlak({1,5,3,7}));
    f.vlakken.emplace_back(Vlak({5,0,6,3}));
    f.vlakken.emplace_back(Vlak({6,2,7,3}));
    f.vlakken.emplace_back(Vlak({0,5,1,4}));

    int curr_point_ind = f.points.size();

    Figuur cylinder = createCylinder(1,10, Color(0,0,0));
    applyTransformation(cylinder,transformationMatrix(0.25, M_PI / 2, 0, 0, Vector3D(0,0,0, true)));
    Figuur cylinder2 = Figuur(cylinder);
    Figuur cylinder3 = Figuur(cylinder);
    Figuur cylinder4 = Figuur(cylinder);
    Figuur cylinder5 = Figuur(cylinder);
    Figuur cylinder6 = Figuur(cylinder);
    applyTransformation(cylinder, translate(Vector3D(1,-1,-1,true)));
    applyTransformation(cylinder2, translate(Vector3D(1,0,-1,true)));
    applyTransformation(cylinder3, translate(Vector3D(1,1,-1,true)));
    applyTransformation(cylinder4, translate(Vector3D(-1,-1,-1,true)));
    applyTransformation(cylinder5, translate(Vector3D(-1,0,-1,true)));
    applyTransformation(cylinder6, translate(Vector3D(-1,1,-1,true)));

    f.addfigure(cylinder);
    f.addfigure(cylinder2);
    f.addfigure(cylinder3);
    f.addfigure(cylinder4);
    f.addfigure(cylinder5);
    f.addfigure(cylinder6);
    return f;
}

Figuur createMengerSpons(Color c, int nr_it){
    Figuur spons = Figuur(c);
    Figuur cube = createCube(c);
    Figures3D temp = {};
    cube.addpoint(0,-1,-1);
    cube.addpoint(0,1,-1);
    cube.addpoint(1,0,-1);
    cube.addpoint(-1,0,-1);
    cube.addpoint(0,-1,1);
    cube.addpoint(0,1,1);
    cube.addpoint(1,0,1);
    cube.addpoint(-1,0,1);
    cube.addpoint(1,-1,0);
    cube.addpoint(-1,-1,0);
    cube.addpoint(-1,1,0);
    cube.addpoint(1,1,0);

    generateFractal(cube,temp,nr_it,3,0);
    for (const auto & fig:temp){
        spons.addfigure(fig);
    }
    return spons;

}


#endif //ENGINE_PREDMADEFIGURES_H
