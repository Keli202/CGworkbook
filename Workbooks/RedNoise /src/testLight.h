
#ifndef REDNOISE_TESTLIGHT_H
#define REDNOISE_TESTLIGHT_H

#endif //REDNOISE_TESTLIGHT_H
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

//bool diffuse=false;
//bool specular=false;
//bool gouraudShading=false;
//bool phongShading=false;

RayTriangleIntersection Get_closest_intersection(vec3 cameraPosition, vec3 rayDirection, const std::vector<ModelTriangle>& triangles) {
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
        if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0 && t > 0) {
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
bool IsShadow(RayTriangleIntersection intersection, vec3 lightPosition,const vector<ModelTriangle>& triangles){
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
            if(t < length(lightPosition-intersection.intersectionPoint)&& i!=intersection.triangleIndex) {
                return true;
            }
        }
    }
    return false;
}

glm::vec3 calculateNormal(const ModelTriangle& triangle,RayTriangleIntersection &intersection) {
    glm::vec3 edge1 = triangle.vertices[1] - triangle.vertices[0];
    glm::vec3 edge2 = triangle.vertices[2] - triangle.vertices[0];
    glm::vec3 normal = glm::cross(edge1, edge2);
    intersection.intersectedTriangle.normal= glm::normalize(normal);
    return intersection.intersectedTriangle.normal;
}

float getProximityLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint) {
    float distance = glm::length(lightPosition - intersectionPoint);
    return std::min(1.0f, std::max(0.0f, 100.0f / static_cast<float>(4 * M_PI * distance * distance)));
}

float getIncidenceLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint, const glm::vec3& normal) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    return std::max(0.0f, glm::dot(glm::normalize(normal), lightDirection));
    //return glm::clamp(glm::dot(glm::normalize(normal), lightDirection), 0.0f, 1.0f);
}

float getSpecularIntensity(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint,  const glm::vec3& normal, float specularExponent,const glm::vec3 cameraPosition) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    glm::vec3 viewDirection = glm::normalize(cameraPosition-intersectionPoint);
    glm::vec3 reflectedRay = (-lightDirection) - (2.0f * normal) * (glm::dot(glm::normalize(-lightDirection), glm::normalize(normal)));
    float reflectionAngle = glm::dot(glm::normalize(viewDirection), glm::normalize(reflectedRay));
    reflectionAngle=glm::max(0.0f,reflectionAngle);
    float SpecularLighting = pow(reflectionAngle, specularExponent);
    return SpecularLighting;
}

float getBrightness(const RayTriangleIntersection &intersection,const glm::vec3 &lightPosition,const glm::vec3 &normal,const glm::vec3 cameraPosition,float specularExponent){
    //proximity
    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);

    //incidence
    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,normal);

    //get brightness
    float brightness = (ProximityLighting * IncidenceLighting);

    //ambient
    brightness = glm::clamp(brightness, 0.2f, 1.0f);

    return brightness;
}

glm::vec3 getRayDirection(float x,float y,glm::vec3 camera_position,mat3 Camera_Orientation){
    float scalingFactor = 150.0f;
    float focal_length = 2.0f;

    double x1 = -(x- camera_position.x - (WIDTH / 2)) / (scalingFactor * focal_length);
    double y1 = -(y-camera_position.y - (HEIGHT / 2)) / (scalingFactor * focal_length);
    double z1 = -focal_length ;

    glm::vec3 ray_direction = glm::vec3(x1 , y1 , z1);
    glm::vec3 direction= normalize(Camera_Orientation*ray_direction);
    return direction;
}

void DrawRaytrace(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation,vec3 &lightPosition,float specularExponent) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            glm::vec3 rayDirection= getRayDirection(x,y,camera_position,Camera_Orientation);
            RayTriangleIntersection intersection = Get_closest_intersection(camera_position, rayDirection, triangles);
            if (!isinf(intersection.distanceFromCamera)) {
                bool inShadow = IsShadow(intersection, lightPosition, triangles);
                Colour colour = intersection.intersectedTriangle.colour;
                glm::vec3 normal = calculateNormal(intersection.intersectedTriangle,intersection);
                float brightness = getBrightness(intersection,lightPosition,normal,camera_position,specularExponent);
                float SpecularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,normal,specularExponent,camera_position);

                if(diffuse){
                    colour.red *=std::min(brightness ,1.0f);
                    colour.green *=std::min(brightness ,1.0f);
                    colour.blue *=std::min(brightness ,1.0f);
                }else if(specular){
                    colour.red = std::min(static_cast<int>(colour.red * brightness+SpecularLighting * 255), 255);
                    colour.green = std::min(static_cast<int>(colour.green * brightness+SpecularLighting * 255), 255);
                    colour.blue = std::min(static_cast<int>(colour.blue * brightness+SpecularLighting * 255), 255);
                }

                if (inShadow) {
                    uint32_t color_value = (255 << 24) + (static_cast<int>(colour.red) / 3 << 16) +
                                           (static_cast<int>(colour.green) / 3 << 8) +
                                           static_cast<int>(colour.blue) / 3;
                    window.setPixelColour(x, y, color_value);
                } else {
                    uint32_t color_value =
                            (255 << 24) + (static_cast<int>(colour.red) << 16) + (static_cast<int>(colour.green) << 8) +
                            static_cast<int>(colour.blue);
                    window.setPixelColour(x, y, color_value);
                }
            }
        }
    }
}


vector<ModelTriangle> getVertexNormals(vector<ModelTriangle> &triangles) {
    for (ModelTriangle& triangle : triangles) {
        glm::vec3 normal = calculateNormal(triangle);
        triangle.normal = glm::normalize(normal);
    }

    for (ModelTriangle& triangle : triangles) {
        for (int i = 0; i < 3; ++i) {
            glm::vec3 vertexNormal = glm::vec3(0.0f);
            int count = 0;
            // Find triangles that share the same vertices
            for (ModelTriangle& sharedTriangle : triangles) {
                    glm::vec3 vertex = triangle.vertices[i];
                    if (vertex == sharedTriangle.vertices[0]||vertex == sharedTriangle.vertices[1]||vertex == sharedTriangle.vertices[2]) {
                        vertexNormal += sharedTriangle.normal;
                        count++;
                }
//                for (int j = 0; j < 3; ++j) {
//                    if (triangle.vertices[i] == sharedTriangle.vertices[j]) {
//                        vertexNormal += sharedTriangle.normal;
//                        count++;
//                        break;
//                    }
//                }
            }
            if (count > 0) {
                glm::vec3 averageNormal = glm::normalize(vertexNormal / static_cast<float>(count));
                triangle.vertexNormals[i] = averageNormal;

            }
        }
    }
    return triangles;
}






glm::vec3 getInterpolatedNormal(const ModelTriangle &triangle, float u, float v) {
    float w = 1.0f - u - v;
    glm::vec3 interpolatedNormal = w * triangle.vertexNormals[0] + u * triangle.vertexNormals[1] + v * triangle.vertexNormals[2];
    return glm::normalize(interpolatedNormal);
}




float getVertexBrightness(const RayTriangleIntersection &intersection, const glm::vec3 &lightPosition, const glm::vec3 &normal, const glm::vec3 &cameraPosition, float specularExponent) {
    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);

    //incidence
    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,normal);

    float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,normal,specularExponent,cameraPosition);

    //get brightness
    float ambientIntensity = 0.1f;


    float brightness = ambientIntensity + ProximityLighting*IncidenceLighting + specularLighting;

    brightness = glm::clamp(brightness, 0.0f, 1.0f);

    return brightness;
}



Colour gouraudShading(const RayTriangleIntersection &intersection, const glm::vec3 &lightPosition, const glm::vec3 &cameraPosition, float specularExponent) {
    const ModelTriangle &triangle = intersection.intersectedTriangle;


    std::vector<float> vertexLightIntensities(3);
    for (int i = 0; i < 3; ++i) {
        glm::vec3 vertex = triangle.vertices[i];
        glm::vec3 vertexNormal = triangle.vertexNormals[i];

        float vertexBrightness = getVertexBrightness(intersection, lightPosition, vertexNormal, cameraPosition, specularExponent);
        vertexLightIntensities[i] = vertexBrightness;
    }

    float u = intersection.u;
    float v = intersection.v;
    float w = 1.0f - u - v;
    float interpolatedLightIntensity = w * vertexLightIntensities[0] + u * vertexLightIntensities[1] + v * vertexLightIntensities[2];

    interpolatedLightIntensity=glm::clamp(interpolatedLightIntensity,0.0f,1.0f);
    Colour colour = intersection.intersectedTriangle.colour;
    colour.red = static_cast<uint8_t>(colour.red * interpolatedLightIntensity);
    colour.green = static_cast<uint8_t>(colour.green * interpolatedLightIntensity);
    colour.blue = static_cast<uint8_t>(colour.blue * interpolatedLightIntensity);
    //cout<<"c:"<<colour<<endl;

    return colour;
}






Colour phongShading(RayTriangleIntersection &intersection, glm::vec3 &lightPosition, glm::vec3 &cameraPosition, float specularExponent) {
    glm::vec3 interpolatedNormal = getInterpolatedNormal(intersection.intersectedTriangle, intersection.u, intersection.v);

    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);
    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,interpolatedNormal);
    float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,interpolatedNormal,specularExponent,cameraPosition);

    Colour colour = intersection.intersectedTriangle.colour;
    colour.red = static_cast<uint8_t>(std::min(colour.red * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
    colour.green = static_cast<uint8_t>(std::min(colour.green * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
    colour.blue = static_cast<uint8_t>(std::min(colour.blue * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));

    return colour;
}


void DrawGouraud(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation,vec3 &lightPosition,float specularExponent) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            glm::vec3 rayDirection= getRayDirection(x,y,camera_position,Camera_Orientation);
            RayTriangleIntersection intersection = Get_closest_intersection(camera_position, rayDirection, triangles);
            if (!isinf(intersection.distanceFromCamera)) {
                bool inShadow = IsShadow(intersection, lightPosition, triangles);
                Colour colour;
                if(gouraud){
                     colour = gouraudShading(intersection, lightPosition, camera_position, specularExponent);
                }else if(phong){
                    colour = phongShading(intersection,lightPosition,camera_position,specularExponent);
                }
                 //cout<<"c:"<<colour<<endl;
                if (inShadow) {
                    uint32_t color_value = (255 << 24) + (static_cast<int>(colour.red) / 3 << 16) +(static_cast<int>(colour.green) / 3 << 8) +static_cast<int>(colour.blue) / 3;
                    window.setPixelColour(x, y, color_value);
                } else {
                    uint32_t color_value = (255 << 24) + (static_cast<int>(colour.red) << 16) + (static_cast<int>(colour.green) << 8) +static_cast<int>(colour.blue);
                    window.setPixelColour(x, y, color_value);
                }
            }
        }
    }
}





