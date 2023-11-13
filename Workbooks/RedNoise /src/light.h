

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


RayTriangleIntersection getClosestIntersection(vec3 cameraPosition, vec3 rayDirection, const std::vector<ModelTriangle>& triangles) {
    RayTriangleIntersection closestIntersection;
    closestIntersection.distanceFromCamera = numeric_limits<float>::infinity();
    //closestIntersection.triangleIndex = -1;
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
                closestIntersection=RayTriangleIntersection(intersectionPoint,closestIntersection.distanceFromCamera,triangle,closestIntersection.triangleIndex);

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
    return std::min(1.0f, std::max(0.0f, 1.0f / static_cast<float>(4 * M_PI * distance * distance)));
}
float calculateIncidenceLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint, const glm::vec3& normal) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    return std::max(0.0f, glm::dot(glm::normalize(normal), lightDirection));
}

float calculateSpecularIntensity(const glm::vec3& lightDirection, const glm::vec3& viewDirection, const glm::vec3& normal, float specularExponent) {
    glm::vec3 reflectDirection = glm::reflect(-lightDirection, normal);
    float specularity = glm::dot(reflectDirection, viewDirection);
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

Colour DiffuseColor(Colour colour,const glm::vec3& lightPosition, const glm::vec3& intersectionPoint, float incidenceLighting){
    float proximityLighting = calculateProximityLighting(lightPosition, intersectionPoint);
    colour.red *= std::min(proximityLighting + incidenceLighting, 1.0f);
    colour.green *= std::min(proximityLighting + incidenceLighting, 1.0f);
    colour.blue *= std::min(proximityLighting + incidenceLighting, 1.0f);

    return colour;
}

Colour SpecularColor(Colour colour,const glm::vec3& lightDirection, const glm::vec3& viewDirection, const glm::vec3& normal, float specularExponent,float incidenceLighting){
    float specularIntensity = calculateSpecularIntensity(lightDirection, viewDirection, normal, specularExponent);
//    colour.red = std::min(255, static_cast<int>(colour.red*incidenceLighting+specularIntensity * 255));
//    colour.green = std::min(255, static_cast<int>(colour.green*incidenceLighting+specularIntensity * 255));
//    colour.blue = std::min(255, static_cast<int>(colour.blue*incidenceLighting+specularIntensity * 255));

    colour.red = std::min(255, static_cast<int>(colour.red+specularIntensity * 255));
    colour.green = std::min(255, static_cast<int>(colour.green+specularIntensity * 255));
    colour.blue = std::min(255, static_cast<int>(colour.blue+specularIntensity * 255));
    const float ambientLightIntensity = 0.2f;
    colour.red = glm::max(static_cast<float>(colour.red), ambientLightIntensity);
    colour.green = glm::max(static_cast<float>(colour.green), ambientLightIntensity);
    colour.blue = glm::max(static_cast<float>(colour.blue), ambientLightIntensity);
    return colour;
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
            //vec3 lightPosition(0.0f, 0.5f, 0.5f);
            if (!isinf(intersection.distanceFromCamera)) {
                bool Shadow = inShadow(intersection, lightPosition, triangles);
                intersection.intersectedTriangle.normal = calculateNormal(intersection.intersectedTriangle);

                float incidenceLighting = calculateIncidenceLighting(lightPosition, intersection.intersectionPoint, intersection.intersectedTriangle.normal);
                glm::vec3 lightDirection = glm::normalize(lightPosition - intersection.intersectionPoint);
                glm::vec3 viewDirection = glm::normalize(camera_position - intersection.intersectionPoint);
                //Colour colour=DiffuseColor(intersection.intersectedTriangle.colour,lightPosition, intersection.intersectionPoint, incidenceLighting);

                //Colour colour=SpecularColor(intersection.intersectedTriangle.colour,lightDirection, viewDirection, intersection.intersectedTriangle.normal, 256.0f);
                Colour colors=DiffuseColor(intersection.intersectedTriangle.colour,lightPosition, intersection.intersectionPoint,incidenceLighting);
                Colour colour=SpecularColor(colors,lightDirection, viewDirection, intersection.intersectedTriangle.normal, 256.0f,incidenceLighting);


                if (Shadow) {
                    uint32_t color_value = (255 << 24) + (static_cast<int>(colour.red) / 3 << 16) +
                                           (static_cast<int>(colour.green) / 3 << 8) +
                                           static_cast<int>(colour.blue) / 3;
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
