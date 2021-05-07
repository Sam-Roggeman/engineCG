#include "src/easy_image.h"
#include "src/ini_configuration.h"
#include "src/vector3d.h"
#include <fstream>
#include <iostream>
#include <string>
#include "list"
#include "src/Line2D.h"
#include "src/Vlak.h"
#include "src/l_parser.h"
#include "src/Figuur.h"
#include "src/helpfunctions.h"
#include <cmath>
#include <limits>
#include <stack>
#include <algorithm>
#include "src/predmadeFigures.h"




using Lines2D = std::list<Line2D>;
using Figures3D = std::list<Figuur>;




img::EasyImage draw2DLines(Lines2D &lines, const int size, const Color &backgroundColor, bool zbuf) {
    double Imagex;
    double Imagey;
    double d;
    double dx;
    double dy;
    calculateScaleOffset(lines,size,d,Imagex,Imagey,dx,dy);
    for (auto &line: lines){
        line.reScale(d);
    }

    for (auto &line: lines){
        line.move(dx, dy);
    }

    auto image = img::EasyImage(roundToInt(Imagex),roundToInt(Imagey));
    ZBuffer zbuffer = ZBuffer(roundToInt(Imagex),roundToInt(Imagey));

    image.clear(backgroundColor.imageColor());
    for (const Line2D &line:lines){
        if (zbuf){
            image.draw_zbuf_line(zbuffer,roundToInt(line.getP1().getX()), roundToInt(line.getP1().getY()),
                                  roundToInt(line.getP2().getX()), roundToInt(line.getP2().getY()),
                                  line.getImageColor(),line.getZ1(),line.getZ2());
        }
        else {
            image.draw_line(roundToInt(line.getP1().getX()), roundToInt(line.getP1().getY()),
                         roundToInt(line.getP2().getX()), roundToInt(line.getP2().getY()),
                         line.getImageColor());
        }
    }
    return image;
}



void toPolar(const Vector3D &point, double &alpha, double & beta, double &r){
      r = sqrt(std::pow(point.x, 2)+std::pow(point.y, 2)+std::pow(point.z, 2));
      alpha = std::atan2(point.y, point.x);
      beta = std::acos(point.z/r);
}



void applyTransformation(Figures3D & figuren, const Matrix& m){
    for (auto &figuur:figuren) {
        applyTransformation(figuur,m);
    }
}





void eyePointTrans(const Vector3D &eyepoint, Figures3D & figuren){
    double alpha;
    double beta;
    double r;
    toPolar(eyepoint,alpha,beta, r);
    Matrix m = eyePointTransformationMatrix(alpha, beta,r);
    applyTransformation(figuren,  m);
}

void eyePointVistrum(const Vector3D &eyepoint, const Vector3D &eye_dir, Figures3D & figuren){
    double alpha;
    double beta;
    double r;
    toPolar(-eye_dir,alpha,beta, r);
    applyTransformation(figuren, translate(-eyepoint)* rotateZ(-(alpha+M_PI/2))* rotateX(-beta));
}


void twoDLSystem(const ini::Configuration &configuration, Lines2D& lines){
    std::string filename = configuration["2DLSystem"]["inputfile"].as_string_or_die();
    Color col = Color(configuration["2DLSystem"]["color"].as_double_tuple_or_die());
    LParser::LSystem2D parseObject;
    std::ifstream input_stream(filename);
    input_stream >> parseObject;
    input_stream.close();
    std::stack<std::tuple<double,double,double>> stack = std::stack<std::tuple<double,double,double>>();

    double angle = parseObject.get_starting_angle();
    std::string finalString = parseObject.get_finalstring();
    double alpha = parseObject.get_angle();
    double currentX = 0;
    double currentY = 0;
    for (char c : finalString ) {
        if (c == '+') {
            angle = (angle + alpha);
            if (angle >= 360) angle -= 360;
        } else if (c == '-') {
            angle = (angle - alpha);
            if (angle < 0) angle += 360;

        } else if (c == '(') stack.push(std::make_tuple(currentX,currentY, angle));
        else if (c == ')') {currentX = std::get<0>(stack.top());
            currentY = std::get<1>(stack.top());
            angle = std::get<2>(stack.top());
            stack.pop();

        }
        else {
            double prevX = currentX;
            double prevY = currentY;
            currentX += cos(angle* M_PI/180);
            currentY += sin(angle* M_PI/180);
            if (parseObject.draw(c)) {
                Point2D p0 = Point2D(prevX,prevY);
                Point2D p1 = Point2D(currentX,currentY);
                lines.push_back(Line2D(p0,p1,col));
            }
        }
    }
}


Point2D doProjection(const Vector3D &point, const double d){
    return {(d*point.x)/(-point.z),(d*point.y)/(-point.z)};
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
                    auto line = Line2D(*point_it_prev, *(point_it_next), figuur.color);
                    lines.emplace_back(line);
                    point_it_prev = point_it_next;
                }
                auto line = Line2D(*points.begin(),*(points.end()-1),figuur.color);
                lines.emplace_back(line);
            }else {
                auto line = Line2D(*points.begin(), *(points.begin()+1), figuur.color);
                lines.emplace_back(line);
            }
            points = {};
        }
    }
    return lines;
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

Figuur threeDLsystem(Color col, const std::string& filename) {
    Figuur f = Figuur(col);


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
        }

        else if (c == '(') {
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
            } else if (c == '|') {
                for (int i = 0; i < 3; i++) {
                    H[i] = -H[i];
                    L[i] = -L[i];
                }
            }
        }
        else {
            std::cerr << "char: '" << c << "' not recognized" << std::endl;
        }
    }
    return f;
}

std::list<Figuur> makeFigures(const ini::Configuration &configuration) {
    Figuur f = Figuur(Color(1,1,1));
    std::list<Figuur> figuren ;
    int nr_figures = configuration["General"]["nrFigures"].as_int_or_die();
    for (int figure_nr = 0; figure_nr < nr_figures; figure_nr++) {
        std::string name = "Figure" + std::to_string(figure_nr);
        std::string type = configuration[name]["type"].as_string_or_die();
        Color c = Color(configuration[name]["color"].as_double_tuple_or_die());

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
                f = createCube(c);
            } else if (type == "Tetrahedron") {
                f = createTetrahedron(c);
            } else if (type == "Octahedron") {
                f = createOctahedron(c);
            } else if (type == "Icosahedron") {
                f = createIcosahedron(c);
            } else if (type == "Dodecahedron") {
                f = createDodecahedron(c);
            } else if (type == "Cylinder") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                f = createCylinder(h, n, c);
            } else if (type == "Cone") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                f = createCone(h, n, c);
            } else if (type == "Sphere") {
                int n = configuration[name]["n"].as_int_or_die();
                f = createSphere(n, c);
            } else if (type == "Torus") {
                double r = configuration[name]["r"].as_double_or_die();
                double R = configuration[name]["R"].as_double_or_die();
                int n = configuration[name]["n"].as_int_or_die();
                int m = configuration[name]["m"].as_int_or_die();
                f = createTorus(r, R, n, m, c);
            } else if (type == "3DLSystem") {
                std::string filename = configuration[name]["inputfile"].as_string_or_die();
                f = threeDLsystem(c, filename);
            }
            else if (type == "BuckyBall") {
                f = createBuckyball(c);
            }
            else if (type == "MengerSponge") {
                int n = configuration[name]["nrIterations"].as_int_or_die();
                f = createMengerSpons(c,n);
            }
                else if (type == "Vrachtwagen"){
                f = createVrachtwagen(c);
            }
            else {
                std::cerr << "type: " << type << " not found" << std::endl;
            }
            applyTransformation(f,
                                transformationMatrix(s, r_x * M_PI / 180, r_y * M_PI / 180, r_z * M_PI / 180, center));
            figuren.emplace_back(f);
        } else {
            if (type == "FractalCube") {
                f = createCube(c);

            } else if (type == "FractalTetrahedron") {
                f = createTetrahedron(c);

            } else if (type == "FractalOctahedron") {
                f = createOctahedron(c);
            } else if (type == "FractalIcosahedron") {
                f = createIcosahedron(c);
            } else if (type == "FractalDodecahedron") {
                f = createDodecahedron(c);
            } else if (type == "FractalCylinder") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                f = createCylinder(h, n, c);
            } else if (type == "FractalCone") {
                int n = configuration[name]["n"].as_int_or_die();
                double h = configuration[name]["height"].as_double_or_die();
                f = createCone(h, n, c);
            } else if (type == "FractalSphere") {
                int n = configuration[name]["n"].as_int_or_die();
                f = createSphere(n, c);
            } else if (type == "FractalTorus") {
                double r = configuration[name]["r"].as_double_or_die();
                double R = configuration[name]["R"].as_double_or_die();
                int n = configuration[name]["n"].as_int_or_die();
                int m = configuration[name]["m"].as_int_or_die();
                f = createTorus(r, R, n, m, c);
            } else if (type == "FractalBuckyBall") {
                f = createBuckyball(c);
            } else if (type == "FractalMengerSponge") {
                int n = configuration[name]["nrIterations"].as_int_or_die();
                f = createMengerSpons(c,n);
            }
            int nr_it = configuration[name]["nrIterations"];
            double fractalScale = configuration[name]["fractalScale"];

            f = generateFractal(f,  nr_it, fractalScale);
            applyTransformation(f,
                                transformationMatrix(s, r_x * M_PI / 180, r_y * M_PI / 180, r_z * M_PI / 180, center));
            figuren.emplace_back(f);
        }
    }
    return figuren;
}

void wireFrame(const ini::Configuration &configuration, Lines2D& lines){
    std::vector<int> eye_pos = configuration["General"]["eye"].as_int_tuple_or_die();
    Vector3D eye_vector = Vector3D(eye_pos[0], eye_pos[1],eye_pos[2], true);
    Figures3D figuren = makeFigures(configuration);

    eyePointTrans(eye_vector,figuren);
    auto lines2 = doProjection(figuren);
    lines.insert(lines.end(), lines2.begin(), lines2.end() );
}

void draw_zbuf_triangle(ZBuffer& zbuf, img::EasyImage& image, Color color,
                        Vector3D const& A, Vector3D const& B, Vector3D const& C,
                        double d, double dx, double dy){

    Point2D proA = Point2D((d*A.x)/(-1*(A.z))+dx, (d*A.y)/(-A.z)+dy);
    Point2D proB = Point2D((d*B.x)/(-1*(B.z))+dx, (d*B.y)/(-B.z)+dy);
    Point2D proC = Point2D((d*C.x)/(-1*(C.z))+dx, (d*C.y)/(-C.z)+dy);
    int min_y = roundToInt(std::min({proA.getY(),proB.getY(),proC.getY()})+0.5);
    int max_y = roundToInt(std::max({proA.getY(),proB.getY(),proC.getY()})-0.5);

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
    Vector3D u = Vector3D(B - A);
    Vector3D v = Vector3D(C - A);
    Vector3D w = v.cross_equals(u);
    double k = w.x* A.x + w.y * A.y + A.z * w.z;
    double dzdx = w.x/(-d*k);
    double dzdy = w.y/(-d*k);


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
void clipPane(std::vector<Vector3D> pane, int i, std::vector<std::vector<Vector3D>>& l) {
    std::vector<std::vector<Vector3D>> nw_l = {};
    double dval;
    double dval2;
    bool a_inside ;
    bool b_inside ;
    bool c_inside ;
    const Vector3D* A;
    const Vector3D* B;
    const Vector3D* C;
    const Vector3D* X;
    const Vector3D* Y;
    Vector3D A2;
    Vector3D B2;
    Vector3D C2;
    Vector3D D2;
    bool near = i == 0;
    bool far = i == 1;
    bool rechts = i ==2;
    bool links = i ==3;
    bool boven = i ==4;
    bool onder = i ==5;
    double dNear;
    double p;
    for ( auto triangle_it = l.begin(); triangle_it < l.end();triangle_it++){
        int count = 0;
        A = &triangle_it->at(0);
        B = &triangle_it->at(1);
        C = &triangle_it->at(2);
        if (near || far) {
            dval = pane[0].z;
            if (near) {
                a_inside = A->z <= dval;
                b_inside = B->z <= dval;
                c_inside = C->z <= dval;
            } else {
                a_inside = A->z >= dval;
                b_inside = B->z >= dval;
                c_inside = C->z >= dval;
            }
            if (a_inside) count++;
            if (b_inside) count++;
            if (c_inside) count++;
            if (count ==3) {
                A2 = *A;
                B2 = *B;
                C2 = *C;
            } else if (count ==0) {
                break;
            } else if (count == 1) {
                if (!a_inside && !b_inside) {
                    //AC
                    p = (dval - C->z) / (A->z - C->z);
                    A2 = (p * *A) + ((1 - p) * *C);
                    //BC
                    p = (dval - C->z) / (B->z - C->z);
                    B2 = (p * *B) + ((1 - p) * *C);
                    C2 = *C;
                } else if (!a_inside && !c_inside) {
                    //AB
                    p = (dval - B->z) / (A->z - B->z);
                    A2 = (p * *A) + ((1 - p) * *B);
                    //CB
                    p = (dval - B->z) / (C->z - B->z);
                    B2 = (p * *C) + ((1 - p) * *B);
                    C2 = *B;
                } else if (!c_inside && !b_inside) {
                    //CA
                    p = (dval - C->z) / (A->z - C->z);
                    A2 = (p * *C) + ((1 - p) * *A);
                    //BA
                    p = (dval - A->z) / (B->z - A->z);
                    B2 = (p * *B) + ((1 - p) * *A);
                    C2 = *A;
                }
            }else if (count == 2){
                if (!a_inside) {
                    //AC
                    X = A; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AB
                    X = A; Y = B;
                    p = (dval - Y->z) / (X->z - Y->z);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *C;
                    C2 = *B;
                } else if (!c_inside) {
                    //AC
                    X = A; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *B;
                } else if (!b_inside) {
                    //AB
                    X = A; Y = B;
                    p = (dval - Y->z) / (X->z - Y->z);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    
                    B2 = *A;
                    C2 = *C;
                }
            }
        }else if (rechts||links){
            dNear = -std::max(std::max(pane[0].z, pane[1].z),std::max(pane[2].z, pane[3].z));
            dval = pane[0].x*(-dNear)/pane[0].z;
            if (rechts) {
                dval2 = A->x*(-dNear)/A->z;
                a_inside = (dval2 <= dval);
                dval2 = B->x*(-dNear)/B->z;
                b_inside = (dval2 <= dval);
                dval2 = C->x*(-dNear)/C->z;
                c_inside = (dval2 <= dval);
            } else {
                dval2 = A->x*(-dNear)/A->z;
                a_inside = (dval2 >= dval);
                dval2 = B->x*(-dNear)/B->z;
                b_inside = (dval2 >= dval);
                dval2 = C->x*(-dNear)/C->z;
                c_inside = (dval2 >= dval);
            }
            if (a_inside) count++;
            if (b_inside) count++;
            if (c_inside) count++;
            if (count ==3) {
                A2 = *A;
                B2 = *B;
                C2 = *C;
            } else if (count ==0) {
                break;
            }else if (count == 1) {
                if (!a_inside && !b_inside) {
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *C;
                } else if (!a_inside && !c_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //CB
                    X = C; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *B;
                } else if (!c_inside && !b_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *A;
                }
            }else if (count == 2){
                if (!a_inside) {
//                    A2 = *A;
//                    B2 = *B;
//                    C2 = *C;
//                    D2 = *C;
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *B;
                    C2 = *C;
                } else if (!c_inside) {
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *B;
                } else if (!b_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *C;
                }
            }
        } else if (boven||onder) {
            dNear = -std::max(std::max(pane[0].z, pane[1].z), std::max(pane[2].z, pane[3].z));
            dval = -pane[0].y * dNear / pane[0].z;
            if (boven) {
                a_inside =  -A->y*dNear/A->z <= dval;
                b_inside = -B->y*dNear/B->z <= dval;
                c_inside = -C->y*dNear/C->z <= dval;
            } else {
                a_inside = -A->y*dNear/A->z >= dval;
                b_inside = -B->y*dNear/B->z >= dval;
                c_inside = -C->y*dNear/C->z >= dval;
            }
            if (a_inside) count++;
            if (b_inside) count++;
            if (c_inside) count++;
            if (count ==3) {
                A2 = *A;
                B2 = *B;
                C2 = *C;
            } else if (count ==0) {
                break;
            } else if (count == 1) {
                if (!a_inside && !b_inside) {
                    //AC
                    X = A;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *C;
                } else if (!a_inside && !c_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //CB
                    X = C;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *B;
                } else if (!c_inside && !b_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *A;
                }
            } else if (count == 2) {
                if (!a_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *B;
                    C2 = *C;
                } else if (!c_inside) {
                    //BC
                    X = B;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //CA
                    X = C;
                    Y = A;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *B;
                    C2 = *A;
                } else if (!b_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *C;
                }
            }
        }
        if (count == 2) {
            nw_l.emplace_back(std::vector<Vector3D>({A2, D2, C2}));
        }
        nw_l.emplace_back(std::vector<Vector3D>({A2,B2,C2}));
    }
    l.clear();
    l = nw_l;
}

std::vector<std::vector<Vector3D>> clipVistrum(const Vector3D &A ,const Vector3D &B,const Vector3D &C,const Figuur &vistrum) {
    std::vector<std::vector<Vector3D>> l;
    l.emplace_back(std::vector<Vector3D>({A,B,C}));
    for (int i = 0; i <6; i++){
        clipPane(vistrum.vlakInPoints(i), i,l);
    }
    return l;
}

img::EasyImage ZBufferTriangles(const ini::Configuration &configuration, const int size) {
    std::vector<int> eye_pos = configuration["General"]["eye"].as_int_tuple_or_die();
    Color bg = Color(configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
    bool clipping = configuration["General"]["clipping"].as_bool_or_default(false);
    std::list<Figuur> figuren = makeFigures(configuration);
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
    eyePointVistrum(eye_vector, eye_dir,figuren);

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
                draw_zbuf_triangle(z, image, figuur.color, figuur.points[vlak.point_indexes[0]], figuur.points[vlak.point_indexes[1]], figuur.points[vlak.point_indexes[2]], d, dx, dy);
            }
        }
    }
    else {
        for (auto &figuur:figuren) {
            for (auto &vlak:figuur.vlakken) {
                triangulated_vakken = triangulate(vlak);
                for (auto &tri_vlak:triangulated_vakken) {
                    draw_zbuf_triangle(z, image, figuur.color, figuur.points[tri_vlak.point_indexes[0]], figuur.points[tri_vlak.point_indexes[1]], figuur.points[tri_vlak.point_indexes[2]], d, dx, dy);
                }
            }
        }
    }

    return image;
}





img::EasyImage generate_image(const ini::Configuration &configuration)
{
    img::EasyImage image;
    std::string type = configuration["General"]["type"].as_string_or_die();
    Color bgcolor = Color(configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
    int size = configuration["General"]["size"].as_int_or_die();
    Lines2D lines = {};
    if (type == "2DLSystem"){
        twoDLSystem(configuration, lines);
        return draw2DLines(lines, size, bgcolor, false);
    }
    else if (type == "Wireframe"){
        wireFrame(configuration,lines);
        return draw2DLines(lines, size, bgcolor, false);
    }
    else if (type == "ZBufferedWireframe"){
        wireFrame(configuration,lines);
        return draw2DLines(lines, size, bgcolor, true);
    }
    else if(type == "ZBuffering"){
        return ZBufferTriangles(configuration, size);
    }
    else {
        std::cerr<< "type: " <<type<< " not recognized"<<std::endl;
    }
}




int main(int argc, char const* argv[])
{
    int retVal = 0;
    try
    {
        for(int i = 1; i < argc; ++i)
        {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(argv[i]);
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                std::string fileName(argv[i]);
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;
                    std::cout << fileName<<" done"<<std::endl;

                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << argv[i] << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}
