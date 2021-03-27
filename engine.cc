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
#include "src/matrices.h"
#include "src/predmadeFigures.h"



using Lines2D = std::list<Line2D>;
using Figures3D = std::list<Figuur>;




img::EasyImage draw2DLines(Lines2D &lines, const int size, const Color &backgroundColor, bool zbuf) {
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
        wireFrame(configuration,lines);
        return draw2DLines(lines, size, bgcolor, true);
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
