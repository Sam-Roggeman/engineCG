//
// Created by User on 8/05/2021.
//

#ifndef ENGINE_ZBUF_H
#define ENGINE_ZBUF_H
#include "3D.h"
#include "Clipping.h"
void draw_zbuf_triangle(ZBuffer& zbuf, img::EasyImage& image,
                        Vector3D const& A, Vector3D const& B, Vector3D const& C,
                        double d, double dx, double dy,
                        const Color ambientReflection, const Color diffuseReflection, const Color specularReflection,
                        const double reflectionCoeff, const Lights3D& lights){

    Point2D proA = Point2D((d*A.x)/(-1*(A.z))+dx, (d*A.y)/(-A.z)+dy);
    Point2D proB = Point2D((d*B.x)/(-1*(B.z))+dx, (d*B.y)/(-B.z)+dy);
    Point2D proC = Point2D((d*C.x)/(-1*(C.z))+dx, (d*C.y)/(-C.z)+dy);
    int min_y = roundToInt(std::min({proA.getY(),proB.getY(),proC.getY()})+0.5);
    int max_y = roundToInt(std::max({proA.getY(),proB.getY(),proC.getY()})-0.5);
    Color color = Color(0,0,0);

    double xlAB;
    double xlAC;
    double xlBC;
    double xrAB;
    double xrAC;
    double xrBC;
    int xr;
    int xl;
    double z_inv;
    double x_g = (proA.getX()+proB.getX()+proC.getX())/3;
    double y_g = (proA.getY()+proB.getY()+proC.getY())/3;
    double z_inv_g = 1 / (3 * A.z) + 1 / (3 * B.z) + 1 / (3 * C.z) ;

    int count;
    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D::cross(u,v);
    double k = w.x* A.x + w.y * A.y + A.z * w.z;
    double dzdx = w.x/(-d*k);
    double dzdy = w.y/(-d*k);
    w.normalise();
    for (const auto & light:lights){
        Vector3D ld;
        if (light->isInfLight()){
            ld = light->getLdVector();
            ld = -ld;
            ld.normalise();
            k =w.dot(ld);
            if (k>0) {
                color += (diffuseReflection * light->diffuseLight * k);
            }
        }
        color +=  (light->ambientLight * ambientReflection);
    }
    for (int y = min_y; y <= max_y; y++){
        xlAB = std::numeric_limits<int>::max();
        xlAC = std::numeric_limits<int>::max();
        xlBC = std::numeric_limits<int>::max();
        xrAB = std::numeric_limits<int>::min();
        xrAC = std::numeric_limits<int>::min();
        xrBC = std::numeric_limits<int>::min();
        count = 0;
        if (proA.getY() != proB.getY() && (y-proA.getY())*(y-proB.getY()) <= 0){
            xrAB = proB.getX() + (proA.getX()-proB.getX())*((y-proB.getY())/(proA.getY()-proB.getY()));
            xlAB = xrAB;
            count +=1;
        }
        if (proA.getY() != proC.getY() && (y-proA.getY())*(y-proC.getY()) <= 0){
            xrAC = proC.getX() + (proA.getX()-proC.getX())*((y-proC.getY())/(proA.getY()-proC.getY()));
            xlAC = xrAC;
            count +=1;
        }
        if (proB.getY() != proC.getY() && (y-proB.getY())*(y-proC.getY()) <= 0){
            xrBC = proC.getX() + (proB.getX()-proC.getX())*((y-proC.getY())/(proB.getY()-proC.getY()));
            xlBC = xrBC;
            count +=1;
        }
        xl = roundToInt(std::min({xlAB,xlBC,xlAC})+0.5);
        xr = roundToInt(std::max({xrAB,xrBC,xrAC})-0.5);
        if (count <2) continue;
        for (int x = xl; x <= xr; x++) {
            z_inv = (z_inv_g) + (x-x_g)*dzdx + (y-y_g)*dzdy;
            if (zbuf.changeIfCloser(y,x,z_inv)){
                image.draw_pixel((unsigned int) x,(unsigned int)y, color.imageColor());
            }
        }
    }
}

img::EasyImage ZBufferTriangles(const ini::Configuration &configuration, const int size) {
    Lights3D lights;
    int nrlights = configuration["General"]["nrLights"].as_int_or_default(0);
    if (nrlights == 0){
        lights.emplace_back(new Light(Color(1,1,1),Color(0,0,0),Color(0,0,0)));
    }
    else{
        std::string light_name;
        for (int light_nr = 0; light_nr < nrlights; light_nr++){
            light_name = "Light"+std::to_string(light_nr);
            if(configuration[light_name]["infinity"].as_bool_or_default(false)){
                std::vector<double> dir = configuration[light_name]["direction"].as_double_tuple_or_die();
                lights.emplace_back(new InfLight(Color(configuration[light_name]["ambientLight"].as_double_tuple_or_default({0,0,0}))
                        ,Color(configuration[light_name]["diffuseLight"].as_double_tuple_or_default({0,0,0}))
                        ,Color(configuration[light_name]["b"].as_double_tuple_or_default({0,0,0}))
                        ,Vector3D::vector(dir[0],dir[1],dir[2])));
            }
            else {
                lights.emplace_back(new Light(Color(configuration[light_name]["ambientLight"].as_double_tuple_or_default({0,0,0}))
                        ,Color(configuration[light_name]["diffuseLight"].as_double_tuple_or_default({0,0,0}))
                        ,Color(configuration[light_name]["b"].as_double_tuple_or_default({0,0,0}))));
            }
        }
    }
    std::vector<int> eye_pos = configuration["General"]["eye"].as_int_tuple_or_die();
    Color bg = Color(configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
    bool clipping = configuration["General"]["clipping"].as_bool_or_default(false);
    std::list<Figuur> figuren;
    makeFigures(configuration,figuren);
    std::list<Figuur> temp_figuren;
    std::vector<Vlak> triangulated_vakken = {};
    double Imagex;
    double Imagey;
    double d;
    double dx;
    double dy;
    std::vector<int> eye_dir_v = configuration["General"]["viewDirection"].as_int_tuple_or_default({-eye_pos[0],-eye_pos[1], -eye_pos[2]});
    Vector3D eye_dir = Vector3D(eye_dir_v[0], eye_dir_v[1],eye_dir_v[2], true);
    Vector3D eye_vector = Vector3D(eye_pos[0], eye_pos[1],eye_pos[2], true);
    eyePointVistrum(eye_vector, eye_dir,figuren,lights);
    Vector3D A;
    Vector3D B;
    Vector3D C;
    std::vector<std::vector<Vector3D>> clipped_triangles;
    Figuur viewVistrum = Figuur(bg);
    if (clipping){
        std::vector<std::vector<Vector3D>> completed;
        viewVistrum = createViewFrustum(configuration["General"]["hfov"].as_double_or_die(),configuration["General"]["dNear"].as_double_or_die(),
                                        configuration["General"]["dFar"].as_double_or_die(), configuration["General"]["aspectRatio"].as_double_or_die(),
                                        bg);
        for (auto &figuur:figuren) {
            completed.clear();
            for (auto &vlak:figuur.vlakken) {
                triangulated_vakken = triangulate(vlak);
                for (auto &tri_vlak:triangulated_vakken) {
                    clipped_triangles = clipVistrum(figuur.points[tri_vlak.point_indexes[0]], figuur.points[tri_vlak.point_indexes[1]], figuur.points[tri_vlak.point_indexes[2]], viewVistrum);
                    completed.reserve(clipped_triangles.size());
                    completed.insert(completed.end(), clipped_triangles.begin(),clipped_triangles.end());
                }
            }
            figuur.points.clear();
            figuur.vlakken.clear();
            for (auto& vlak: completed){
                figuur.addVlak(vlak);
            }
        }
    }
    calculateScaleOffset(doProjection(figuren),size,d,Imagex,Imagey,dx,dy);
    img::EasyImage image = img::EasyImage(roundToInt(Imagex), roundToInt(Imagey));
    image.clear(bg.imageColor());
    ZBuffer z = ZBuffer(roundToInt(Imagex), roundToInt(Imagey));
    if (clipping){
        for (const auto &figuur:figuren) {
            for (const auto &vlak:figuur.vlakken) {
                draw_zbuf_triangle(z, image,figuur.points[vlak.point_indexes[0]], figuur.points[vlak.point_indexes[1]], figuur.points[vlak.point_indexes[2]], d, dx, dy,figuur.ambientReflection,figuur.diffuseReflection,figuur.specularReflection,figuur.reflectionCoefficient,lights);
            }
        }
    }
    else {
        for (auto &figuur:figuren) {
            for (auto &vlak:figuur.vlakken) {
                triangulated_vakken = triangulate(vlak);
                for (auto &tri_vlak:triangulated_vakken) {
                    draw_zbuf_triangle(z, image, figuur.points[tri_vlak.point_indexes[0]], figuur.points[tri_vlak.point_indexes[1]], figuur.points[tri_vlak.point_indexes[2]], d, dx, dy,figuur.ambientReflection,figuur.diffuseReflection,figuur.specularReflection,figuur.reflectionCoefficient,lights);
                }
            }
        }
    }
    for (auto& light:lights){
        delete light;
    }
    return image;
}
#endif //ENGINE_ZBUF_H
