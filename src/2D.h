//
// Created by User on 8/05/2021.
//

#ifndef ENGINE_2D_H
#define ENGINE_2D_H
#include "helpfunctions.h"
#include "stack"
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
#endif //ENGINE_2D_H
