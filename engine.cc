#include "easy_image.h"
#include "ini_configuration.h"
#include "vector3d.h"
#include <fstream>
#include <iostream>
#include <string>
#include "list"
#include "Line2D.h"
#include "Vlak.h"
#include "l_parser.h"
#include "Figuur.h"
#include <cmath>
#include <limits>
#include <stack>



using Lines2D = std::list<Line2D>;
using Figures3D = std::list<Figuur>;

inline int roundToInt(double d)
{
    return static_cast<int>(std::round(d));
}


void findExtrema(double &xmin, double &xmax, double &ymin, double &ymax, const Lines2D& lines){
    if ( lines.empty()){std::cerr<<"No lines" << std::endl;}
    for (const auto &line:lines){
        Point2D p = line.getP1();
        for (int i = 0; i<2;i++){
            if (p.getX() < xmin) xmin = p.getX();
            else if (p.getX() > xmax) xmax = p.getX();
            if (p.getY() < ymin) ymin = p.getY();
            else if (p.getY() > ymax) ymax = p.getY();
            p = line.getP2();
        }
    }
}


img::EasyImage *draw2DLines(Lines2D &lines, const int size, const Color &backgroundColor, bool zbuf) {
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();
    findExtrema(xmin, xmax, ymin, ymax, lines);
    double xrange = xmax-xmin;
    double yrange = ymax-ymin;
    double Imagex = size*((xrange)/(fmax(xrange,yrange)));
    double Imagey = size*((yrange)/(fmax(xrange,yrange)));
    double d = 0.95*(Imagex/xrange);
    for (auto &line: lines){
        line.reScale(d);
    }
    double DCx = d*((xmin+xmax)/2);
    double DCy = d*((ymin+ymax)/2);
    double dx = Imagex/2 - DCx;
    double dy = Imagey/2 - DCy;
    for (auto &line: lines){
        line.move(dx, dy);
    }
    auto* image = new img::EasyImage(roundToInt(Imagex),roundToInt(Imagey));
    ZBuffer zbuffer = ZBuffer(roundToInt(Imagex),roundToInt(Imagey));

    image->clear(backgroundColor.imageColor());
    for (const Line2D &line:lines){
        if (zbuf){
            image->draw_zbuf_line(zbuffer,roundToInt(line.getP1().getX()), roundToInt(line.getP1().getY()),
                                  roundToInt(line.getP2().getX()), roundToInt(line.getP2().getY()),
                                  line.getImageColor(),line.getZ1(),line.getZ2());
        }
        else {
            image->draw_line(roundToInt(line.getP1().getX()), roundToInt(line.getP1().getY()),
                         roundToInt(line.getP2().getX()), roundToInt(line.getP2().getY()),
                         line.getImageColor());
        }
    }
    return image;
}

Matrix scalingMatrix(const double scale){
    Matrix S = Matrix();
    for (int i = 1; i <= 3; i++)
        S(i,i) = scale;
    return S;
}

Matrix rotateX(const double alpha){
    Matrix m_x = Matrix();
    m_x(2,2) = cos(alpha);
    m_x(3,3) = cos(alpha);
    m_x(3,2) = -sin(alpha);
    m_x(2,3) = sin(alpha);
    return m_x;
}

Matrix rotateY(const double alpha){
    Matrix m_y = Matrix();
    m_y(1,1) = cos(alpha);
    m_y(1,3) = -sin(alpha);
    m_y(3,1) = sin(alpha);
    m_y(3,3) = cos(alpha);
    return m_y;
}

Matrix rotateZ(const double alpha){
    Matrix m_z = Matrix();
    m_z(1,1) = cos(alpha);
    m_z(2,2) = cos(alpha);
    m_z(2,1) = -sin(alpha);
    m_z(1,2) = sin(alpha);
    return m_z;
}

Matrix translate(const Vector3D &vector){
    Matrix m_t = Matrix();
    m_t(4,1)= vector.x;
    m_t(4,2)= vector.y;
    m_t(4,3)= vector.z;
    return m_t;
}

void toPolar(const Vector3D &point, double &alpha, double & beta, double &r){
      r = sqrt(std::pow(point.x, 2)+std::pow(point.y, 2)+std::pow(point.z, 2));
      alpha = std::atan2(point.y, point.x);
      beta = std::acos(point.z/r);
}

Matrix transformationMatrix(const double scale, const double alpha_x, const double alpha_y,const double alpha_z, const Vector3D &vector){
    return scalingMatrix(scale)*rotateX(alpha_x)*rotateY(alpha_y)*rotateZ(alpha_z)*translate(vector);
}

void applyTransformation(Figuur & f, const Matrix & m){
    for (auto &punt : f.points){
        punt *= m;
    }
}

void applyTransformation(Figures3D & figuren, const Matrix& m){
    for (auto &figuur:figuren) {
        applyTransformation(figuur,m);
    }
}

Matrix eyePointTransformationMatrix(const double alpha, const double beta, const double r){
    Matrix m = Matrix();
    m(1,1) = -sin(alpha);
    m(1,2) = -cos(alpha)*cos(beta);
    m(1,3) = cos(alpha)*sin(beta);
    m(2,1) = cos(alpha);
    m(2,2) = -sin(alpha)*cos(beta);
    m(2,3) = sin(alpha)*sin(beta);
    m(3,2) = sin(beta);
    m(3,3) = cos(beta);
    m(4,3) = -r;
    return m;
}



void eyePointTrans(const Vector3D &eyepoint, Figures3D & figuren){
    double alpha;
    double beta;
    double r;
    toPolar(eyepoint,alpha,beta, r);
    Matrix m = eyePointTransformationMatrix(alpha, beta,r);
    applyTransformation(figuren,  m);

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
                Point2D p= doProjection(figuur.points[point_index-1], 1.00);
                p.setZ(figuur.points[point_index-1].z);
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
            }
            else {
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
            v.point_indexes.emplace_back(punt_ind+1);
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
    f.vlakken.emplace_back(Vlak({1,5,3,7}));
    f.vlakken.emplace_back(Vlak({5,2,8,3}));
    f.vlakken.emplace_back(Vlak({2,6,4,8}));
    f.vlakken.emplace_back(Vlak({6,1,7,4}));
    f.vlakken.emplace_back(Vlak({7,3,8,4}));
    f.vlakken.emplace_back(Vlak({1,6,2,5}));
    return f;
}

Figuur createTetrahedron(Color color) {
    Figuur f = Figuur(color);
    f.points.emplace_back(Vector3D::point(1,-1,-1));
    f.points.emplace_back(Vector3D::point(-1,1,-1));
    f.points.emplace_back(Vector3D::point(1,1,1));
    f.points.emplace_back(Vector3D::point(-1,-1,1));
    f.vlakken.emplace_back(Vlak({1,2,3}));
    f.vlakken.emplace_back(Vlak({2,4,3}));
    f.vlakken.emplace_back(Vlak({1,4,2}));
    f.vlakken.emplace_back(Vlak({1,3,4}));
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
    f.vlakken.emplace_back(Vlak({1,2,6}));
    f.vlakken.emplace_back(Vlak({2,3,6}));
    f.vlakken.emplace_back(Vlak({3,4,6}));
    f.vlakken.emplace_back(Vlak({4,1,6}));
    f.vlakken.emplace_back(Vlak({2,1,5}));
    f.vlakken.emplace_back(Vlak({3,2,5}));
    f.vlakken.emplace_back(Vlak({4,3,5}));
    f.vlakken.emplace_back(Vlak({1,4,5}));
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

    f.vlakken.emplace_back(Vlak({1,2,3}));
    f.vlakken.emplace_back(Vlak({1,3,4}));
    f.vlakken.emplace_back(Vlak({1,4,5}));
    f.vlakken.emplace_back(Vlak({1,5,6}));
    f.vlakken.emplace_back(Vlak({1,6,2}));

    f.vlakken.emplace_back(Vlak({2,7,3}));
    f.vlakken.emplace_back(Vlak({3,7,8}));
    f.vlakken.emplace_back(Vlak({3,8,4}));
    f.vlakken.emplace_back(Vlak({4,8,9}));
    f.vlakken.emplace_back(Vlak({4,9,5}));

    f.vlakken.emplace_back(Vlak({5,9,10}));
    f.vlakken.emplace_back(Vlak({5,10,6}));
    f.vlakken.emplace_back(Vlak({6,10,11}));
    f.vlakken.emplace_back(Vlak({6,11,2}));
    f.vlakken.emplace_back(Vlak({2,11,7}));

    f.vlakken.emplace_back(Vlak({12,8,7}));
    f.vlakken.emplace_back(Vlak({12,9,8}));
    f.vlakken.emplace_back(Vlak({12,10,9}));
    f.vlakken.emplace_back(Vlak({12,11,10}));
    f.vlakken.emplace_back(Vlak({12,7,11}));

    return f;
}

Figuur createDodecahedron(Color color) {
    Figuur f = Figuur(color);
    Figuur ico = createIcosahedron(color);
    for (const Vlak& ico_vlak:ico.vlakken){
        double x = 0;
        double y = 0;
        double z = 0;
        for (int ico_ind:ico_vlak.point_indexes) {
            x += ico.points[ico_ind-1].x;
            y += ico.points[ico_ind-1].y;
            z += ico.points[ico_ind-1].z;
        }
        f.points.emplace_back(Vector3D::point(x/3,y/3,z/3));
    }
    f.vlakken.emplace_back(Vlak({1,2,3,4,5}));
    f.vlakken.emplace_back(Vlak({1,6,7,8,2}));
    f.vlakken.emplace_back(Vlak({2,8,9,10,3}));
    f.vlakken.emplace_back(Vlak({3,10,11,12,4}));
    f.vlakken.emplace_back(Vlak({4,12,13,14,5}));

    f.vlakken.emplace_back(Vlak({5,14,15,6,1}));
    f.vlakken.emplace_back(Vlak({20,19,18,17,16}));
    f.vlakken.emplace_back(Vlak({20,15,14,13,19}));
    f.vlakken.emplace_back(Vlak({19,13,12,11,18}));
    f.vlakken.emplace_back(Vlak({18,11,10,9,17}));

    f.vlakken.emplace_back(Vlak({17,9,8,7,16}));
    f.vlakken.emplace_back(Vlak({16,7,6,15,20}));

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
            Vector3D p1 = f.points[p_ind1-1];
            Vector3D p2 = f.points[p_ind2-1];
            Vector3D p3 = f.points[p_ind3-1];
            Vector3D p4 = Vector3D::point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2);
            Vector3D p5 = Vector3D::point((p1.x + p3.x) / 2, (p1.y + p3.y) / 2, (p1.z + p3.z) / 2);
            Vector3D p6 = Vector3D::point((p2.x + p3.x) / 2, (p2.y + p3.y) / 2, (p2.z + p3.z) / 2);
            int p_ind4 = f.points.size()+1;
            int p_ind5 = f.points.size()+2;
            int p_ind6 = f.points.size()+3;
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
        f.vlakken.emplace_back(Vlak({i+1, (i+1)%(nr_vlakken)+1,nr_vlakken+1}));
        grd_ind.emplace_back(i+1);
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
        f.vlakken.emplace_back(Vlak({i+1, (i+1)%(nr_vlakken)+1, nr_vlakken+(i+1)%(nr_vlakken)+1,i+1+nr_vlakken}));
        grd_ind.emplace_back(i+1);
        bov_ind.emplace_back(i+1+nr_vlakken);
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
            f.vlakken.emplace_back(Vlak({i*n+j+1, n*((i+1)%n) + j + 1, n*((i+1)%n) + (j+1)%m +1 ,n*i+1+(j+1)%m}));
        }
    }
    return f;
}

Figuur threeDLsystem(Color col, std::string filename) {
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
                f.vlakken.emplace_back(Vlak({i+1,i+2}));
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


void wireFrame(const ini::Configuration &configuration, Lines2D& lines){
    Figuur f = Figuur(Color(1,1,1));
    int nr_figures = configuration["General"]["nrFigures"].as_int_or_die();
    std::vector<int> eye_pos = configuration["General"]["eye"].as_int_tuple_or_die();
    Vector3D eye_vector = Vector3D(eye_pos[0], eye_pos[1],eye_pos[2], true);
    Figures3D figuren = {};
    for (int figure_nr = 0; figure_nr < nr_figures; figure_nr++){
        std::string name = "Figure"+std::to_string(figure_nr);
        std::string type = configuration[name]["type"].as_string_or_die();
        Color c = Color(configuration[name]["color"].as_double_tuple_or_die());

        double s = configuration[name]["scale"].as_double_or_die();
        double r_x = configuration[name]["rotateX"].as_double_or_die();
        double r_y = configuration[name]["rotateY"].as_double_or_die();
        double r_z = configuration[name]["rotateZ"].as_double_or_die();
        std::vector<double> v = configuration[name]["center"].as_double_tuple_or_die();
        Vector3D center = Vector3D::vector(v[0],v[1],v[2]);
        if (type == "LineDrawing"){
            f = lineDrawing(configuration, name);
        }
        else if (type == "Cube"){
            f = createCube(c);
        }
        else if (type == "Tetrahedron"){
            f = createTetrahedron(c);
        }
        else if (type == "Octahedron"){
            f = createOctahedron(c);
        }
        else if (type == "Icosahedron"){
            f = createIcosahedron(c);
        }
        else if (type == "Dodecahedron"){
            f = createDodecahedron(c);
        }
        else if (type == "Cylinder"){
            int n = configuration[name]["n"].as_int_or_die();
            double h = configuration[name]["height"].as_double_or_die();
            f = createCylinder(h,n,c);
        }
        else if (type == "Cone"){
            int n = configuration[name]["n"].as_int_or_die();
            double h = configuration[name]["height"].as_double_or_die();
            f = createCone(h,n,c);
        }
        else if (type == "Sphere"){
            int n = configuration[name]["n"].as_int_or_die();
            f = createSphere(n,c);
        }
        else if (type == "Torus"){
            double r = configuration[name]["r"].as_double_or_die();
            double R = configuration[name]["R"].as_double_or_die();
            int n = configuration[name]["n"].as_int_or_die();
            int m= configuration[name]["m"].as_int_or_die();
            f = createTorus(r, R, n,m, c);
        }
        else if (type == "3DLSystem"){
            std::string filename = configuration[name]["inputfile"].as_string_or_die();
            f = threeDLsystem(c, filename);
        }
        else {
            std::cerr << "type: " << type << " not found"<<std::endl;
        }
        applyTransformation(f,transformationMatrix(s,r_x*M_PI/180,r_y*M_PI/180,r_z*M_PI/180,center));
        figuren.emplace_back(f);

    }
    eyePointTrans(eye_vector,figuren);
    auto lines2 = doProjection(figuren);
    lines.insert(lines.end(), lines2.begin(), lines2.end() );
}



img::EasyImage generate_image(const ini::Configuration &configuration)
{
    std::string type = configuration["General"]["type"].as_string_or_die();
    Color bgcolor = Color(configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
    int size = configuration["General"]["size"].as_int_or_die();
    Lines2D lines = {};
    if (type == "2DLSystem"){
        twoDLSystem(configuration, lines);
        return *draw2DLines(lines, size, bgcolor, false);
    }
    else if (type == "Wireframe"){
        wireFrame(configuration,lines);
        return *draw2DLines(lines, size, bgcolor, false);
    }
    else if (type == "ZBufferedWireframe"){
        wireFrame(configuration,lines);
        return *draw2DLines(lines, size, bgcolor, true);
    }
    else if(type == "ZBuffering"){
        wireFrame(configuration,lines);
        return *draw2DLines(lines, size, bgcolor, true);
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
