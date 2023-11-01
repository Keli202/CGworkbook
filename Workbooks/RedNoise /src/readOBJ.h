#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <ModelTriangle.h>

using namespace std;
using namespace glm;

#ifndef REDNOISE_READOBJ_H
#define REDNOISE_READOBJ_H

#endif //REDNOISE_READOBJ_H
void readOBJColour(const std::string& filename,std::vector<std::pair<std::string, Colour>>& palette) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open material file " << filename << std::endl;
        return;
    }
    string line;
    //vector<vec3> vertices;
    string content;
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0]=="newmtl"){
            //content +=line+ '\n';
            //palette.push_back(std::make_pair(currentMaterial, currentColour));
            content = tokens[1];
        }
        if(tokens[0]=="Kd"){
            float r, g, b;
            r = std::stof(tokens[1]) ;
            g = std::stof(tokens[2]) ;
            b = std::stof(tokens[3]) ;
            int red = static_cast<int>(r * 255);
            int green = static_cast<int>(g * 255);
            int blue = static_cast<int>(b * 255);
            palette.emplace_back(content, Colour(red, green, blue));
        }

    }
    file.close();
}

vector<ModelTriangle> readOBJ(const std::string& filename, float scale){
    ifstream file(filename);
    vector<ModelTriangle> triangles;
    vector<vec3> vertices;

    if(!file.is_open()){
        cerr<<"Error: Could not open file "<<filename<<std::endl;
        return triangles;
    }

    string line;

    // std::vector<Colour> colours;
    vector<pair<string, Colour>> palette;

    // Read the MTL file and populate the palette
    readOBJColour("../cornell-box.mtl",palette);
    string currentMaterial;

    while(getline(file,line)){
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "v") {
            vec3 vertex;
            //cout<<tokens.size()<<endl;
            vertex.x= std::stof(tokens[1]);
            vertex.y = std::stof(tokens[2]);
            vertex.z = std::stof(tokens[3]);
            vertex *= scale;
            vertices.push_back(vertex);
        }
        if (tokens[0] == "f") {
            std::vector<std::string> parts = split(line, '/');
            int v1, v2, v3;
            v1 = std::stoi(tokens[1]) - 1;
            v2 = std::stoi(tokens[2]) - 1;
            v3 = std::stoi(tokens[3]) - 1;

            Colour triangleColour; // Default colour
            for (const auto& data : palette) {
                if (data.first == currentMaterial) {
                    triangleColour = data.second;
                    //cout<<currentMaterial<<endl;
                    break;
                }
            }
            //cout<<triangleColour<<endl;
            ModelTriangle triangle(vertices[v1], vertices[v2], vertices[v3], triangleColour);
            triangles.push_back(triangle);
//               ModelTriangle triangle(vertices[v1], vertices[v2], vertices[v3],Colour{255,255,255});
//               triangles.push_back(triangle);
        }
        if(tokens[0] == "usemtl"){
            currentMaterial = tokens[1];
        }
    }

    for (const ModelTriangle& triangle : triangles) {
        std::cout << triangle << std::endl;
    }

    file.close();
    return triangles;
//return vertices;
}