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








void readOBJColour(const std::string& filename, std::vector<std::pair<std::string, Colour>>& palette) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open material file " << filename << std::endl;
        return;
    }
    std::string line;
    std::string currentMaterial;
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "newmtl") {
            currentMaterial = tokens[1];
        } else if (tokens[0] == "Kd") {
            float r = std::stof(tokens[1]);
            float g = std::stof(tokens[2]);
            float b = std::stof(tokens[3]);
            int red = static_cast<int>(r * 255);
            int green = static_cast<int>(g * 255);
            int blue = static_cast<int>(b * 255);
            palette.emplace_back(currentMaterial, Colour(red, green, blue));
        } else if (tokens[0] == "map_Kd") {

            palette.back().second.name = tokens[1];
        }
    }
    file.close();
}
std::vector<ModelTriangle> readOBJ(const std::string& filename, float scale) {
    std::ifstream file(filename);
    std::vector<ModelTriangle> triangles;
    std::vector<vec3> vertices;
    std::vector<TexturePoint> texturePoints;
    std::vector<std::pair<std::string, Colour>> palette;

    readOBJColour("../textured-cornell-box.mtl", palette);

    std::string currentMaterial;
    std::string line;

    while (getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "v") {
            vec3 vertex(std::stof(tokens[1]) * scale, std::stof(tokens[2]) * scale, std::stof(tokens[3]) * scale);
            vertices.push_back(vertex);
        } else if (tokens[0] == "vt") {
            TexturePoint texPoint(std::stof(tokens[1]), std::stof(tokens[2]));
            texturePoints.push_back(texPoint);
        } else if (tokens[0] == "f") {
            std::vector<std::string> v1 = split(tokens[1], '/');
            std::vector<std::string> v2 = split(tokens[2], '/');
            std::vector<std::string> v3 = split(tokens[3], '/');
            Colour triangleColour(255, 255, 255);
            for (const auto& item : palette) {
                if (item.first == currentMaterial) {
                    triangleColour = item.second;
                    break;
                }
            }
            ModelTriangle triangle(
                    vertices[stoi(v1[0]) - 1],
                    vertices[stoi(v2[0]) - 1],
                    vertices[stoi(v3[0]) - 1],
                    triangleColour
            );
            if (!v1[1].empty() && !v2[1].empty() && !v3[1].empty() && v1.size() > 1 && v2.size() > 1 && v3.size() > 1) {
                triangle.texturePoints[0] = texturePoints[stoi(v1[1]) - 1];
                triangle.texturePoints[1] = texturePoints[stoi(v2[1]) - 1];
                triangle.texturePoints[2] = texturePoints[stoi(v3[1]) - 1];
            }
            //triangle.normal= calculateNormal(triangle);
            triangles.push_back(triangle);
        } else if (tokens[0] == "usemtl") {
            currentMaterial = tokens[1];
        }
    }

    file.close();
    return triangles;
}
