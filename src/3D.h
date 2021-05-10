//
// Created by User on 8/05/2021.
//

#ifndef ENGINE_3D_H
#define ENGINE_3D_H
#include "stack"
#include "predmadeFigures.h"
#include "helpfunctions.h"
#include "l_parser.h"


Point2D doProjection(const Vector3D &point, const double &d){
    return {(d*point.x)/(-point.z),(d*point.y)/(-point.z)};
}
Figuur lineDrawing(const ini::Configuration &configuration, const std::string& figure_name) {
    Figuur fig = Figuur(Color(configuration[figure_name]["color"].as_double_tuple_or_die()));
    int nr_points = configuration[figure_name]["nrPoints"].as_int_or_die();
    int nr_lines = configuration[figure_name]["nrLines"].as_int_or_die();

    for (int line_nr = 0; line_nr < nr_lines; line_nr++) {
        std::string line_name = "line" + std::to_string(line_nr);
        Vlak v = Vlak();
        auto punt_inds = configuration[figure_name][line_name].as_int_tuple_or_die();
        for (auto const &punt_ind:punt_inds) {
            v.point_indexes.emplace_back(punt_ind);
        }
        fig.vlakken.emplace_back(v);
    }
    for (int point_nr = 0; point_nr < nr_points; point_nr++) {
        std::string point_name = "point" + std::to_string(point_nr);
        auto punt_pos = configuration[figure_name][point_name].as_double_tuple_or_die();
        Vector3D v = Vector3D();
        v.x = punt_pos[0];
        v.y = punt_pos[1];
        v.z = punt_pos[2];
        fig.points.emplace_back(v);
    }
    return fig;
}
Figuur threeDLsystem(Figuur& f, const std::string& filename) {
    LParser::LSystem3D parseObject;
    std::ifstream input_stream(filename);
    input_stream >> parseObject;
    input_stream.close();
    std::stack<std::tuple<std::vector<double>,std::vector<double>,std::vector<double>, double,double,double>> stack =
            std::stack<std::tuple<std::vector<double>,std::vector<double>,std::vector<double>, double,double,double>>();
    std::vector<double> H = {1,0,0};
    std::vector<double> U = {0,0,1};
    std::vector<double> L = {0,1,0};
    std::vector<double> H_temp_copy;
    std::vector<double> L_temp_copy;
    std::vector<double> U_temp_copy;
    const double alpha = parseObject.get_angle();
    std::string finalString = parseObject.get_finalstring();
    double x = 0;
    double y = 0;
    double z = 0;
    double prevx;
    double prevy;
    double prevz;
    int count = 0;
    for (char c : finalString ) {
        count +=1;
        if (parseObject.cInAlphabet(c)){
            prevx = x;
            prevy = y;
            prevz = z;
            x+= H[0];
            y+= H[1];
            z+= H[2];
            if (parseObject.draw(c)) {
                int i = f.points.size();
                f.addpoint(prevx, prevy, prevz);
                f.addpoint(x, y, z);
                f.vlakken.emplace_back(Vlak({i,i+1}));
            }
        }else if (c == '(') {
            stack.push(std::make_tuple(H,L,U,x,y,z));
        }
        else if (c == ')') {
            H = std::get<0>(stack.top());
            L = std::get<1>(stack.top());
            U = std::get<2>(stack.top());
            x = std::get<3>(stack.top());
            y = std::get<4>(stack.top());
            z = std::get<5>(stack.top());
            stack.pop();
        }
        else if (c == '+' || c == '-'|| c=='^'||c =='&'||c=='\\'||c=='/'|| c =='|') {
            H_temp_copy = H;
            L_temp_copy = L;
            U_temp_copy = U;
            if (c == '+') {
                for (int i = 0; i < 3; i++) {
                    H[i] = H_temp_copy[i] * cos(alpha*(M_PI/180)) + L_temp_copy[i] * sin(alpha*(M_PI/180));
                    L[i] = -H_temp_copy[i] * sin(alpha*(M_PI/180)) + L_temp_copy[i] * cos(alpha*(M_PI/180));
                }
            } else if (c == '-') {
                for (int i = 0; i < 3; i++) {
                    H[i] = H_temp_copy[i] * cos(-alpha*(M_PI/180)) + L_temp_copy[i] * sin(-alpha*(M_PI/180));
                    L[i] = -H_temp_copy[i] * sin(-alpha*(M_PI/180)) + L_temp_copy[i] * cos(-alpha*(M_PI/180));
                }
            } else if (c == '^') {
                for (int i = 0; i < 3; i++) {
                    H[i] = H_temp_copy[i] * cos(alpha*(M_PI/180)) + U_temp_copy[i] * sin(alpha*(M_PI/180));
                    U[i] = -H_temp_copy[i] * sin(alpha*(M_PI/180)) + U_temp_copy[i] * cos(alpha*(M_PI/180));
                }
            } else if (c == '&') {
                for (int i = 0; i < 3; i++) {
                    H[i] = H_temp_copy[i] * cos(-alpha*(M_PI/180)) + U_temp_copy[i] * sin(-alpha*(M_PI/180));
                    U[i] = -H_temp_copy[i] * sin(-alpha*(M_PI/180)) + U_temp_copy[i] * cos(-alpha*(M_PI/180));
                }
            } else if (c == '\\') {
                for (int i = 0; i < 3; i++) {
                    L[i] = (L_temp_copy[i] * cos(alpha*(M_PI/180)) - U_temp_copy[i] * sin(alpha*(M_PI/180)));
                    U[i] = (L_temp_copy[i] * sin(alpha*(M_PI/180)) + U_temp_copy[i] * cos(alpha*(M_PI/180)));
                }
            } else if (c == '/') {
                for (int i = 0; i < 3; i++) {
                    L[i] = L_temp_copy[i] * cos(-alpha*(M_PI/180)) - U_temp_copy[i] * sin(-alpha*(M_PI/180));
                    U[i] = L_temp_copy[i] * sin(-alpha*(M_PI/180)) + U_temp_copy[i] * cos(-alpha*(M_PI/180));
                }
            } else {
                for (int i = 0; i < 3; i++) {
                    H[i] = -H[i];
                    L[i] = -L[i];
                }
            }
        }else {
            std::cerr << "char: '" << c << "' not recognized" << std::endl;
        }
    }
    return f;
}

void makeFigures(const ini::Configuration &configuration, Figures3D& figuren) {
    Figuur f = Figuur();
    Color am;
    Color diff = Color(0,0,0);
    Color spec = Color(0,0,0);
    double refflco = 0;
    int nr_figures = configuration["General"]["nrFigures"].as_int_or_die();
    for (int figure_nr = 0; figure_nr < nr_figures; figure_nr++) {
        std::string name = "Figure" + std::to_string(figure_nr);
        std::string type = configuration[name]["type"].as_string_or_die();
        if (configuration[name]["color"].exists()){
            am = Color(configuration[name]["color"].as_double_tuple_or_die());
        }else{
            am = Color(configuration[name]["ambientReflection"].as_double_tuple_or_default({0,0,0}));
            diff = Color(configuration[name]["diffuseReflection"].as_double_tuple_or_default({0,0,0}));
            spec = Color(configuration[name]["specularReflection"].as_double_tuple_or_default({0,0,0}));
            refflco = configuration[name]["reflectionCoefficient"].as_double_or_default(0);
        }
        f = Figuur(am,diff,spec,refflco);

        double s = configuration[name]["scale"].as_double_or_die();
        double r_x = configuration[name]["rotateX"].as_double_or_die();
        double r_y = configuration[name]["rotateY"].as_double_or_die();
        double r_z = configuration[name]["rotateZ"].as_double_or_die();
        std::vector<double> v = configuration[name]["center"].as_double_tuple_or_die();
        Vector3D center = Vector3D::vector(v[0], v[1], v[2]);
        if (type.size()<8 ||type.substr(0, 7) != "Fractal") {
            if (type == "LineDrawing") {
                f = lineDrawing(configuration, name);
            } else if (type == "Cube") {
                createCube(f);
            } else if (type == "Tetrahedron") {
                createTetrahedron(f);
            } else if (type == "Octahedron") {
                createOctahedron(f);
            } else if (type == "Icosahedron") {
                createIcosahedron(f);
            } else if (type == "Dodecahedron") {
                createDodecahedron(f);
            } else if (type == "Cylinder") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                createCylinder(h, n, f);
            }else if (type == "Cone") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                createCone(h, n, f);
            }else if (type == "Sphere") {
                int n = configuration[name]["n"].as_int_or_die();
                createSphere(n, f);
            }else if (type == "Torus") {
                double r = configuration[name]["r"].as_double_or_die();
                double R = configuration[name]["R"].as_double_or_die();
                int n = configuration[name]["n"].as_int_or_die();
                int m = configuration[name]["m"].as_int_or_die();
                createTorus(r, R, n, m, f);
            }else if (type == "3DLSystem") {
                std::string filename = configuration[name]["inputfile"].as_string_or_die();
                threeDLsystem(f,filename);
            }else if (type == "BuckyBall") {
                createBuckyball(f);
            }else if (type == "MengerSponge") {
                int n = configuration[name]["nrIterations"].as_int_or_die();
                createMengerSpons(f,n);
            }
            else if (type == "Vrachtwagen"){
                createVrachtwagen(f);
            }
            else {
                std::cerr << "type: " << type << " not found" << std::endl;
            }
            applyTransformation(f,
                                transformationMatrix(s, r_x * M_PI / 180, r_y * M_PI / 180, r_z * M_PI / 180, center));
            figuren.emplace_back(f);
        } else {
            if (type == "FractalCube") {
                createCube(f);
            } else if (type == "FractalTetrahedron") {
                createTetrahedron(f);
            } else if (type == "FractalOctahedron") {
                createOctahedron(f);
            } else if (type == "FractalIcosahedron") {
                createIcosahedron(f);
            } else if (type == "FractalDodecahedron") {
                createDodecahedron(f);
            } else if (type == "FractalCylinder") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                createCylinder(h, n, f);
            } else if (type == "FractalCone") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                createCone(h, n, f);
            } else if (type == "FractalSphere") {
                int n = configuration[name]["n"].as_int_or_die();
                createSphere(n, f);
            } else if (type == "FractalTorus") {
                double r = configuration[name]["r"].as_double_or_die();
                double R = configuration[name]["R"].as_double_or_die();
                int n = configuration[name]["n"].as_int_or_die();
                int m = configuration[name]["m"].as_int_or_die();
                createTorus(r, R, n, m, f);
            } else if (type == "FractalBuckyBall") {
                createBuckyball(f);
            } else if (type == "FractalMengerSponge") {
                int n = configuration[name]["nrIterations"].as_int_or_die();
                createMengerSpons(f,n);
            }
            int nr_it = configuration[name]["nrIterations"];
            double fractalScale = configuration[name]["fractalScale"];

            generateFractal(f,  nr_it, fractalScale);
            applyTransformation(f,
                                transformationMatrix(s, r_x * M_PI / 180, r_y * M_PI / 180, r_z * M_PI / 180, center));
            figuren.emplace_back(f);
        }
    }
}

Lines2D doProjection(const Figures3D & figuren){
    std::vector<Point2D> points;
    Lines2D lines;
    for (const auto &figuur:figuren){
        for (const auto& vlak:figuur.vlakken){
            for (int point_index:vlak.point_indexes) {
                Point2D p= doProjection(figuur.points[point_index], 1.00);
                p.setZ(figuur.points[point_index].z);
                points.emplace_back(p);
            }
            if (points.size() > 2) {
                auto point_it_prev = points.begin();
                for (auto point_it_next = points.begin()+1; point_it_next != points.end(); point_it_next++) {
                    auto line = Line2D(*point_it_prev, *(point_it_next), figuur.ambientReflection);
                    lines.emplace_back(line);
                    point_it_prev = point_it_next;
                }
                auto line = Line2D(*points.begin(),*(points.end()-1),figuur.ambientReflection);
                lines.emplace_back(line);
            }else {
                auto line = Line2D(*points.begin(), *(points.begin()+1), figuur.ambientReflection);
                lines.emplace_back(line);
            }
            points = {};
        }
    }
    return lines;
}
void wireFrame(const ini::Configuration &configuration, Lines2D& lines){
    std::vector<int> eye_pos = configuration["General"]["eye"].as_int_tuple_or_die();
    Vector3D eye_vector = Vector3D(eye_pos[0], eye_pos[1],eye_pos[2], true);
    Figures3D figuren;
    makeFigures(configuration,figuren);
    Lights3D l = {};
    eyePointTrans(eye_vector,figuren,l);
    auto lines2 = doProjection(figuren);
    lines.insert(lines.end(), lines2.begin(), lines2.end() );
}


#endif //ENGINE_3D_H
