#include "easy_image.h"
#include "ini_configuration.h"
#include "vector3d.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "list"
#include "Line2D.h"
#include "l_parser.h"
#include "Figuur.h"
#include <math.h>
#include <limits>
#include <stack>

using Lines2D = std::list<Line2D>;


inline int roundToInt(double d)
{
    return static_cast<int>(std::round(d));
}


void findExtrema(double &xmin, double &xmax, double &ymin, double &ymax, const Lines2D& lines){
    if ( lines.empty()){std::cerr<<"No lines";}
    for (Line2D line:lines){
        Point2D* p = line.getP1();
        for (int i = 0; i<2;i++){
            if (p->getX() < xmin) xmin = p->getX();
            if (p->getY() < ymin) ymin = p->getY();
            if (p->getX() > xmax) xmax = p->getX();
            if (p->getY() > ymax) ymax = p->getY();
            p = line.getP2();
        }
    }
}


img::EasyImage* draw2DLines(const Lines2D &lines, const int size, const Color& backgroundColor) {
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
    for (auto line: lines){
        line.reScale(d);
    }
    double DCx = d*((xmin+xmax)/2);
    double DCy = d*((ymin+ymax)/2);
    double dx = Imagex/2 - DCx;
    double dy = Imagey/2 - DCy;
    for (auto line: lines){
        line.move(dx, dy);
    }
    img::EasyImage* image = new img::EasyImage(roundToInt(Imagex),roundToInt(Imagey));
    image->clear(backgroundColor.imageColor());
    for (Line2D line:lines){
        image->draw_line(roundToInt(line.getP1()->getX()), roundToInt(line.getP1()->getY()),
                         roundToInt(line.getP2()->getX()), roundToInt(line.getP2()->getY()),
                         line.getImageColor());
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
    for (auto punt : f.points){
        punt *= m;
    }
}

void applyTransformation(std::vector<Figuur> & figuren, const Matrix& m){
    for (auto figuur:figuren) {
        for (auto punt : figuur.points) {
            punt *= m;
        }
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

void eyePointTrans(const Vector3D &eyepoint, std::vector<Figuur>& figuren){
    double alpha;
    double beta;
    double r;
    toPolar(eyepoint,alpha,beta, r);
    Matrix m = eyePointTransformationMatrix(alpha, beta,r);
    applyTransformation(figuren,  m);

}

void twoDLSystem(const ini::Configuration &configuration, Lines2D& lines){
    std::string filename = configuration["2DLSystem"]["inputfile"].as_string_or_die();
    Color* col = new Color(configuration["2DLSystem"]["color"].as_double_tuple_or_die());
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
                Point2D* p0 = new Point2D(prevX,prevY);
                Point2D* p1 = new Point2D(currentX,currentY);
                lines.push_back(Line2D(p0,p1,col));
            }
        }
    }
}





img::EasyImage generate_image(const ini::Configuration &configuration)
{
    std::string type = configuration["General"]["type"].as_string_or_die();
    Color bgcolor = Color(configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
    int size = configuration["General"]["size"].as_int_or_die();
    Lines2D lines = {};
    if (type == "2DLSystem"){
        twoDLSystem(configuration, lines);
        return *draw2DLines(lines, size, bgcolor);
    }
    else if (type == "Wireframe"){
        Matrix s = scalingMatrix(3);
        s.print(std::cout);
    }
    return img::EasyImage();
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
