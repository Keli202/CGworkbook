

#ifndef REDNOISE_LIGHT_H
#define REDNOISE_LIGHT_H

#endif //REDNOISE_LIGHT_H
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#define WIDTH 320*2.5
#define HEIGHT 240*2.5
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
using namespace std;
using namespace glm;

bool diffuse=false;
bool specular=false;
bool gouraud=false;
bool phong=false;

RayTriangleIntersection getClosestIntersection(vec3 cameraPosition, vec3 rayDirection, const std::vector<ModelTriangle>& triangles) {
    RayTriangleIntersection closestIntersection;
    closestIntersection.distanceFromCamera = numeric_limits<float>::infinity();
    for(int i = 0; i < triangles.size(); i++) {
        ModelTriangle triangle = triangles[i];

        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 spVector = cameraPosition - triangle.vertices[0];
        mat3 DEMatrix(-rayDirection, e0, e1);
        vec3 possibleSolution = inverse(DEMatrix) * spVector;

        float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
        // Check if there is an intersection within the triangle
        if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0 && t > 0) {
            // Check if this intersection is closer than previous ones
            if(t < closestIntersection.distanceFromCamera) {
                vec3 intersectionPoint=cameraPosition+possibleSolution.x*rayDirection;
                closestIntersection.distanceFromCamera = t;
                closestIntersection.triangleIndex = i;
                closestIntersection=RayTriangleIntersection(intersectionPoint,closestIntersection.distanceFromCamera,triangle,closestIntersection.triangleIndex,u,v);

            }
        }
    }

    return closestIntersection;
}
bool inShadow(RayTriangleIntersection intersection, vec3 lightPosition,const vector<ModelTriangle>& triangles){
    vec3 shadow_ray= normalize(lightPosition-intersection.intersectionPoint);
    for(int i = 0; i < triangles.size(); i++) {
        ModelTriangle triangle = triangles[i];

        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 spVector = intersection.intersectionPoint - triangle.vertices[0];
        mat3 DEMatrix(-shadow_ray, e0, e1);
        vec3 possibleSolution = inverse(DEMatrix) * spVector;

        float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
        if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0 &&t>0.01) {
            if(t < length(lightPosition-intersection.intersectionPoint)&&i!=intersection.triangleIndex) {
                return true;
            }
        }
    }
    return false;
}


float calculateProximityLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint) {
    float distance = glm::length(lightPosition - intersectionPoint);
    return std::min(1.0f, std::max(0.0f, 10.0f / static_cast<float>(4 * M_PI * distance * distance)));
}

float calculateIncidenceLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint, const glm::vec3& normal) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    return std::max(0.0f, glm::dot(glm::normalize(normal), lightDirection));
    //return glm::clamp(glm::dot(glm::normalize(normal), lightDirection), 0.0f, 1.0f);
}

float calculateSpecularIntensity(const glm::vec3& lightDirection, const glm::vec3& viewDirection, const glm::vec3& normal, float specularExponent) {
    glm::vec3 reflectDirection = glm::reflect(-lightDirection, normal);
    //glm::vec3 reflectDirection= -lightDirection - 2.0f * normal * glm::dot(glm::normalize(-lightDirection), glm::normalize(normal));
    float specularity = glm::dot(glm::normalize(reflectDirection), glm::normalize(viewDirection));
    //float specularity = glm::dot(glm::normalize(viewDirection), glm::normalize(reflectDirection));

    specularity=glm::max(0.0f,specularity);
    float specularIntensity = pow(specularity, specularExponent);
    return specularIntensity;
}

glm::vec3 calculateNormal(const ModelTriangle& triangle) {
    glm::vec3 edge1 = triangle.vertices[1] - triangle.vertices[0];
    glm::vec3 edge2 = triangle.vertices[2] - triangle.vertices[0];
    glm::vec3 normal = glm::cross(edge1, edge2);
    return glm::normalize(normal);
}



float GetBrightness(const glm::vec3& intersectionPoint,vec3 lightPosition,glm::vec3 normal){
    float incidenceLighting = calculateIncidenceLighting(lightPosition, intersectionPoint, normal);
    float proximityLighting = calculateProximityLighting(lightPosition, intersectionPoint);
    float brightness=glm::clamp((proximityLighting * incidenceLighting), 0.0f, 1.0f);
    return brightness;
}

void drawRaytrace(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation,vec3 lightPosition) {


    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            float scalingFactor = 150.0f;
            float focal_length = 2.0f;

            double x1 = -(x- camera_position.x - (WIDTH / 2)) / (scalingFactor * focal_length);
            double y1 = -(y-camera_position.y - (HEIGHT / 2)) / (scalingFactor * focal_length);
            double z1 = -focal_length ;

            glm::vec3 ray_direction = glm::vec3(x1 , y1 , z1);
            glm::vec3 direction= normalize(Camera_Orientation*ray_direction);
            // Find the closest intersection
            RayTriangleIntersection intersection = getClosestIntersection(camera_position, direction, triangles);
            if (!isinf(intersection.distanceFromCamera)) {
                bool Shadow = inShadow(intersection, lightPosition, triangles);
                intersection.intersectedTriangle.normal = calculateNormal(intersection.intersectedTriangle);
                glm::vec3 normal= intersection.intersectedTriangle.normal;
                glm::vec3 lightDirection = glm::normalize(lightPosition - intersection.intersectionPoint);
                glm::vec3 viewDirection = glm::normalize(camera_position - intersection.intersectionPoint);

                float specularIntensity = calculateSpecularIntensity(lightDirection, viewDirection, intersection.intersectedTriangle.normal, 256.0f);
                float brightness = GetBrightness(intersection.intersectionPoint,lightPosition,intersection.intersectedTriangle.normal);
                brightness = glm::clamp(brightness, 0.20f, 1.0f);
                Colour colour=intersection.intersectedTriangle.colour;
                if(diffuse){
                    //colour=DiffuseColor(intersection.intersectedTriangle.colour,lightPosition, intersection.intersectionPoint,incidenceLighting);
                    colour.red *= std::min(brightness ,1.0f);
                    colour.green *=std::min( brightness ,1.0f);
                    colour.blue *= std::min( brightness ,1.0f);

                }else if(specular){
                    //Highlights are more prominent and can be clearly observed, suitable for very smooth or metallic objects
                    if (specularIntensity > brightness){
                        brightness = specularIntensity;
                    }
                    colour.red = std::min(static_cast<int>(colour.red * brightness), 255);
                    colour.green = std::min(static_cast<int>(colour.green * brightness), 255);
                    colour.blue = std::min(static_cast<int>(colour.blue * brightness), 255);


                    //It more realistically simulates the actual behavior of light on an object's surface
//                    colour.red = std::min(static_cast<int>(colour.red * brightness+specularIntensity * 255), 255);
//                    colour.green = std::min(static_cast<int>(colour.green * brightness+specularIntensity * 255), 255);
//                    colour.blue = std::min(static_cast<int>(colour.blue * brightness+specularIntensity * 255), 255);
                }

                if (Shadow) {
//                    brightness=0.4;
//                    colour.red *= brightness ;
//                    colour.green *= brightness ;
//                    colour.blue *=  brightness ;
                    uint32_t color_value = (255 << 24) + (static_cast<int>(colour.red)/3  << 16) +
                                           (static_cast<int>(colour.green)/3  << 8) +
                                           static_cast<int>(colour.blue)/3  ;
                    window.setPixelColour(x, y, color_value);
                } else{
                    uint32_t color_value =
                            (255 << 24) + (static_cast<int>(colour.red) << 16) + (static_cast<int>(colour.green) << 8) +
                            static_cast<int>(colour.blue);
                    window.setPixelColour(x, y, color_value);
                }
            }
        }
    }
}









//vector<ModelTriangle> calculateVertexNormals(vector<ModelTriangle> &triangles) {
//    for (auto& triangle : triangles) {
//        triangle.normal = calculateFaceNormal(triangle);
//    }
//    for (auto& triangle : triangles) {
//        triangle.vertexNormals.resize(triangle.vertices.size());
//        for (int v = 0; v < triangle.vertices.size(); ++v) {
//            vec3 vertexNormal = vec3(0.0f);
//            int count = 0;
//            for (auto& otherTriangle : triangles) {
//                for (auto& otherVertex : otherTriangle.vertices) {
//                    if (triangle.vertices[v] == otherVertex) {
//                        vertexNormal += otherTriangle.normal;
//                        count++;
//                    }
//                }
//            }
//            if (count != 0) {
//                triangle.vertexNormals[v] = normalize(vertexNormal / static_cast<float>(count));
//            }
//        }
//    }
//    return triangles;
//}
//
//Colour calculateLightingAtVertex(RayTriangleIntersection r, vec3 &normal, vec3 &light, vec3 &vertex, vec3 &cameraPosition, float shininess ) {
//    float ambientStrength = 0.2f;
//    Colour ambientColour = {
//            int(r.intersectedTriangle.colour.red * ambientStrength),
//            int(r.intersectedTriangle.colour.green * ambientStrength),
//            int(r.intersectedTriangle.colour.blue * ambientStrength)
//    };
//
//
//    vec3 toLight = light - r.intersectionPoint;
//    float distance = length(toLight);
//    float Intensity =  100 /(4*M_PI*(pow(distance,2)));
//    float attenuation = std::max(0.0f, std::min(1.0f,Intensity));
//
//    vec3 lightDirection = normalize(toLight);
//    float diff = std::max(dot(normal, lightDirection), 0.0f);
//    Colour diffuseColour = {
//            int(r.intersectedTriangle.colour.red * diff * attenuation),
//            int(r.intersectedTriangle.colour.green * diff * attenuation),
//            int(r.intersectedTriangle.colour.blue * diff * attenuation)
//    };
//
//    vec3 viewDir = normalize(cameraPosition - r.intersectionPoint);
//    vec3 reflectDir = normalize(normalize(-lightDirection) - (normal * 2.0f * dot(normalize(-lightDirection), normal)));
//    float spec = pow(std::max(dot(viewDir, reflectDir), 0.0f), shininess);
//    Colour specularColour = {
//            int(255.0f * spec),
//            int(255.0f *  spec),
//            int(255.0f * spec)
//    };
//
//    int finalRed = std::min(ambientColour.red + diffuseColour.red + specularColour.red, 255);
//    int finalGreen = std::min(ambientColour.green + diffuseColour.green + specularColour.green, 255);
//    int finalBlue = std::min(ambientColour.blue + diffuseColour.blue + specularColour.blue, 255);
//
//    return Colour{finalRed, finalGreen, finalBlue};
//}
//
//
//Colour interpolateVertexColours(ModelTriangle &triangle, vec3 &barycentric) {
//    Colour colour;
//    colour.red = barycentric.x * triangle.vertexColours[0].red + barycentric.y * triangle.vertexColours[1].red + barycentric.z * triangle.vertexColours[2].red;
//    colour.green = barycentric.x * triangle.vertexColours[0].green + barycentric.y * triangle.vertexColours[1].green + barycentric.z * triangle.vertexColours[2].green;
//    colour.blue = barycentric.x * triangle.vertexColours[0].blue + barycentric.y * triangle.vertexColours[1].blue + barycentric.z * triangle.vertexColours[2].blue;
//    return colour;
//}
//
//uint32_t convertColourToUint(Colour &colour) {
//    uint32_t r = static_cast<uint32_t>(std::min(std::max(colour.red, 0), 255));
//    uint32_t g = static_cast<uint32_t>(std::min(std::max(colour.green, 0), 255));
//    uint32_t b = static_cast<uint32_t>(std::min(std::max(colour.blue, 0), 255));
//    return (255 << 24) + (r << 16) + (g << 8) + b;
//}
//
//vec3 getBarycentricCoordinates(vec3 &point, ModelTriangle &triangle) {
//    vec3 v0 = triangle.vertices[1] - triangle.vertices[0];
//    vec3 v1 = triangle.vertices[2] - triangle.vertices[0];
//    vec3 v2 = point - triangle.vertices[0];
//
//    float d00 = dot(v0, v0);
//    float d01 = dot(v0, v1);
//    float d11 = dot(v1, v1);
//    float d20 = dot(v2, v0);
//    float d21 = dot(v2, v1);
//    float denom = d00 * d11 - d01 * d01;
//
//    float v = (d11 * d20 - d01 * d21) / denom;
//    float w = (d00 * d21 - d01 * d20) / denom;
//    float u = 1.0f - v - w;
//
//    return vec3(u, v, w);
//}
//
//void drawGouraud(vector<ModelTriangle> triangles, DrawingWindow &window, vec3 &cameraPosition, float focalLength, mat3 &cameraOrientation,vec3 &light,float scaling_factor){
//    for (int y = 0; y < window.height; y++) {
//        for (int x = 0; x < window.width; x++) {
//            vec3 rayDirection = getRayDirectionFromCanvas(x,y,window.width,window.height,focalLength,scaling_factor,cameraOrientation,cameraPosition);
//            RayTriangleIntersection r = getClosestValidIntersection(triangles, cameraPosition, rayDirection);
//            for (auto& triangle : triangles) {
//                triangle.vertexColours.resize(triangle.vertices.size());
//                for (int v = 0; v < triangle.vertices.size(); ++v) {
//                    triangle.vertexColours[v] = calculateLightingAtVertex(r,triangle.vertexNormals[v],light,triangle.vertices[v],cameraPosition,64.0f);
//                }
//            }
//
//            if (!isinf(r.distanceFromCamera)) {
//                vec3 barycentric = getBarycentricCoordinates(r.intersectionPoint, r.intersectedTriangle);
//                Colour interpolatedColour = interpolateVertexColours(r.intersectedTriangle, barycentric);
//                if (is_shadow(r, triangles, light)) {
//                    interpolatedColour = Colour{int(interpolatedColour.red * 0.5), int(interpolatedColour.green * 0.5), int(interpolatedColour.blue * 0.5)};
//                }
//                window.setPixelColour(x, y, convertColourToUint(interpolatedColour));
//            }
//        }
//    }
//}
//





