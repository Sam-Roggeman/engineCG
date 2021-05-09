#include "src/easy_image.h"
#include "src/ini_configuration.h"
#include "src/vector3d.h"
#include <fstream>
#include <iostream>
#include <string>
#include "list"
#include "src/Line2D.h"
#include "src/Vlak.h"
#include "src/Figuur.h"
#include "src/helpfunctions.h"
#include <algorithm>
#include "src/ZBUF.h"
#include "src/3D.h"
#include "src/2D.h"

using Lines2D = std::list<Line2D>;
using Figures3D = std::list<Figuur>;







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
    else if(type == "LightedZBuffering"){
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
