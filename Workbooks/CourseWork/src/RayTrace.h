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
#ifndef REDNOISE_RAYTRACE_H
#define REDNOISE_RAYTRACE_H
#endif //REDNOISE_RAYTRACE_H
//vector <vec3> lightSources() {
//    vec3 lightPosition(0.0f, 0.6f, 1.0f);
//    int numSamples = 10;
//    float sampleRadius = 0.1f;
//
//    for (int i = 0; i < numSamples; i++) {
//        float angle = (2 * M_PI / numSamples) * i;
//        lightSources().push_back(lightPosition + glm::vec3(sampleRadius * cos(angle), 0, sampleRadius * sin(angle)));
//    }
//    return lightSources();
//}
//vector <vec3> lightSource=lightSources();
vector <vec3> lightSource{

        vec3(0.2f, 0.6f, 1.0f),
        vec3(-0.2f, 0.6f, 1.0f),
        vec3(0.2, 0.6, 1.2),
         vec3(-0.2, 0.6, 1.2),
         vec3(0.0, 0.6, 1.0),
         vec3(0.0, 0.6, 1.2),
         vec3(0.0, 0.6, 0.8),
         vec3(-0.2, 0.6, 0.8),
         vec3(-0.2, 0.6, 1.2),
         //vec3(0.1, 0.6, 0.0),
//         vec3(0.1, 0.8, 0.1),
//         vec3(-0.1, 0.9, -0.1),
//         vec3(-0.1, 0.9, 0.0),
//         vec3(-0.1, 0.9, 0.1),
//         vec3(0.0, 0.9, -0.1),
//         vec3(0.0, 0.9, 0.0),
//         vec3(0.0, 0.9, 0.1),
//         vec3(0.1, 0.9, -0.1),
//         vec3(0.1, 0.9, 0.0),
//         vec3(0.1, 0.9, 0.1),
//        vec3(-0.2, 1.0, 1.8),
//        vec3(-0.2, 1.0, 2.0),
//        vec3(-0.2, 1.0, 2.2),
//        vec3(0.0, 1.0, 1.8),
//        vec3(0.0, 1.0, 2.0),
//        vec3(0.0, 1.0, 2.2),
//        vec3(0.2, 1.0, 1.8),
//        vec3(0.2, 1.0, 2.0),
//        vec3(0.2, 1.0, 2.2)

};










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
bool IsShadow(RayTriangleIntersection intersection, vector<vec3>& lightSources,const vector<ModelTriangle>& triangles){
    for (const auto& lightPosition : lightSources) {
        vec3 shadow_ray = normalize(lightPosition - intersection.intersectionPoint);
        for (int i = 0; i < triangles.size(); i++) {
            ModelTriangle triangle = triangles[i];

            vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
            vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
            vec3 spVector = intersection.intersectionPoint - triangle.vertices[0];
            mat3 DEMatrix(-shadow_ray, e0, e1);
            vec3 possibleSolution = inverse(DEMatrix) * spVector;

            float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
            if (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0 && t > 0.01) {
                if (t < length(lightPosition - intersection.intersectionPoint) && i != intersection.triangleIndex) {
                    return true;
                }
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
glm::vec3 calculateFaceNormal(const ModelTriangle& triangle) {
    glm::vec3 edge1 = triangle.vertices[1] - triangle.vertices[0];
    glm::vec3 edge2 = triangle.vertices[2] - triangle.vertices[0];
    glm::vec3 normal = glm::cross(edge1, edge2);
    return glm::normalize(normal);
}


float getProximityLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint) {
    float distance = glm::length(lightPosition - intersectionPoint);
    return std::min(1.0f, std::max(0.0f, 10.0f / static_cast<float>(4 * M_PI * distance * distance)));
}

float getIncidenceLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint, const glm::vec3& normal) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    return std::max(0.0f, glm::dot(glm::normalize(normal), lightDirection));
    //return glm::clamp(glm::dot(glm::normalize(normal), lightDirection), 0.0f, 1.0f);
}

float getSpecularIntensity(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint,  const glm::vec3& normal, float specularExponent,const glm::vec3 cameraPosition,mat3 Camera_Orientation) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    glm::vec3 viewDirection = glm::normalize(cameraPosition-intersectionPoint);
    //viewDirection=glm::normalize(Camera_Orientation*viewDirection);
    glm::vec3 reflectedRay = glm::reflect(-lightDirection, normal);
    //glm::vec3 reflectedRay = (-lightDirection) - (2.0f * normal) * (glm::dot(glm::normalize(-lightDirection), glm::normalize(normal)));
    // float reflectionAngle = glm::dot(glm::normalize(-lightDirection), glm::normalize(reflectedRay));
    float reflectionAngle = glm::dot(glm::normalize(viewDirection), glm::normalize(reflectedRay));
    reflectionAngle=glm::max(0.0f,reflectionAngle);
    float SpecularLighting = pow(reflectionAngle, specularExponent);
    return SpecularLighting;
}

float getBrightness(const RayTriangleIntersection &intersection,vector<vec3>& lightSources,const glm::vec3 &normal,const glm::vec3 cameraPosition,float specularExponent,const vector<ModelTriangle>& triangles){
//    //proximity
//    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);
//
//    //incidence
//    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,normal);
//
//    //get brightness
//    float brightness = (ProximityLighting * IncidenceLighting);
//
//    //ambient
//    brightness = glm::clamp(brightness/lightSources.size(), 0.2f, 1.0f);
//
//    return brightness;
//}
    float totalBrightness = 0.0f;
    for (const auto& lightPosition : lightSources) {
        if (!IsShadow(intersection, lightSources, triangles)) {
            float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
            float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, normal);
            float brightness = ProximityLighting * IncidenceLighting;
            brightness=glm::clamp(brightness, 0.2f, 1.0f);
            totalBrightness += brightness;
            cout<<"b:"<<totalBrightness<<endl;
        }
    }
    return glm::clamp(totalBrightness / lightSources.size(), 0.2f, 1.0f);
}




vector<ModelTriangle> getVertexNormals(vector<ModelTriangle> &triangles) {
    for ( ModelTriangle& triangle : triangles) {
        glm::vec3 normal = calculateFaceNormal(triangle);
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




float getVertexBrightness(const RayTriangleIntersection &intersection, vector<vec3>& lightSources, const glm::vec3 &normal, const glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation,const vector<ModelTriangle>& triangles) {
//    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);
//
//    //incidence
//    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,normal);
//
//    float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,normal,specularExponent,cameraPosition,Camera_Orientation);
//
//    //get brightness
//    float ambientIntensity = 0.1f;
//
//
//    float brightness = ambientIntensity + ProximityLighting*IncidenceLighting + specularLighting;
//
//    brightness = glm::clamp(brightness/lightSources.size(), 0.0f, 1.0f);


    float totalBrightness = 0.0f;
    for (const auto& lightPosition : lightSources) {
        if (!IsShadow(intersection, lightSources, triangles)) {
            float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
            float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, normal);
            float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,normal,specularExponent,cameraPosition,Camera_Orientation);
            float ambientIntensity = 0.1f;
            float brightness = ambientIntensity + ProximityLighting*IncidenceLighting + specularLighting;
            brightness=glm::clamp(brightness, 0.2f, 1.0f);
            totalBrightness += brightness;
            //cout<<"bv:"<<totalBrightness<<endl;

        }
    }
    return glm::clamp(totalBrightness / lightSources.size(), 0.2f, 1.0f);

    //return brightness;
}



//Colour gouraudShading(const RayTriangleIntersection &intersection, const glm::vec3 &lightPosition, const glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation,const vector<ModelTriangle>& triangles) {
//    const ModelTriangle &triangle = intersection.intersectedTriangle;
//
//
//    std::vector<float> vertexLightIntensities(3);
//    for (int i = 0; i < 3; ++i) {
//        glm::vec3 vertex = triangle.vertices[i];
//        glm::vec3 vertexNormal = triangle.vertexNormals[i];
//
//        float vertexBrightness = getVertexBrightness(intersection, lightPosition, vertexNormal, cameraPosition, specularExponent,Camera_Orientation,triangles);
//        vertexLightIntensities[i] = vertexBrightness;
//    }
//
//    float u = intersection.u;
//    float v = intersection.v;
//    float w = 1.0f - u - v;
//    float interpolatedLightIntensity = w * vertexLightIntensities[0] + u * vertexLightIntensities[1] + v * vertexLightIntensities[2];
//
//    interpolatedLightIntensity=glm::clamp(interpolatedLightIntensity,0.0f,1.0f);
//    Colour colour = intersection.intersectedTriangle.colour;
//    colour.red = static_cast<uint8_t>(colour.red * interpolatedLightIntensity);
//    colour.green = static_cast<uint8_t>(colour.green * interpolatedLightIntensity);
//    colour.blue = static_cast<uint8_t>(colour.blue * interpolatedLightIntensity);
//    //cout<<"c:"<<colour<<endl;
//
//    return colour;
//}


Colour gouraudShading(const RayTriangleIntersection &intersection, vector<vec3>& lightSources, const glm::vec3 &cameraPosition, float specularExponent, mat3 Camera_Orientation, const vector<ModelTriangle>& triangles) {
    const ModelTriangle &triangle = intersection.intersectedTriangle;
    std::vector<float> vertexLightIntensities(3, 0.0f);

    for (const auto& lightPosition : lightSources) {
        for (int i = 0; i < 3; ++i) {
            glm::vec3 vertex = triangle.vertices[i];
            glm::vec3 vertexNormal = triangle.vertexNormals[i];

            float vertexBrightness = getVertexBrightness(intersection, lightSources, vertexNormal, cameraPosition, specularExponent, Camera_Orientation, triangles);
            vertexLightIntensities[i] += vertexBrightness;
        }
    }

    for (float &intensity : vertexLightIntensities) {
        intensity /= lightSources.size();
    }

    float u = intersection.u;
    float v = intersection.v;
    float w = 1.0f - u - v;
    float interpolatedLightIntensity = w * vertexLightIntensities[0] + u * vertexLightIntensities[1] + v * vertexLightIntensities[2];

    interpolatedLightIntensity = glm::clamp(interpolatedLightIntensity, 0.0f, 1.0f);
    Colour colour = intersection.intersectedTriangle.colour;
    colour.red = static_cast<uint8_t>(colour.red * interpolatedLightIntensity);
    colour.green = static_cast<uint8_t>(colour.green * interpolatedLightIntensity);
    colour.blue = static_cast<uint8_t>(colour.blue * interpolatedLightIntensity);

    return colour;
}




//Colour phongShading(RayTriangleIntersection &intersection, glm::vec3 &lightPosition, glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation) {
//    glm::vec3 interpolatedNormal = getInterpolatedNormal(intersection.intersectedTriangle, intersection.u, intersection.v);
//
//    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);
//    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,interpolatedNormal);
//    float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,interpolatedNormal,specularExponent,cameraPosition,Camera_Orientation);
//
//    Colour colour = intersection.intersectedTriangle.colour;
//    colour.red = static_cast<uint8_t>(std::min(colour.red * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
//    colour.green = static_cast<uint8_t>(std::min(colour.green * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
//    colour.blue = static_cast<uint8_t>(std::min(colour.blue * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
//
//    return colour;
//}

Colour phongShading(RayTriangleIntersection &intersection, vector<vec3>& lightSources, glm::vec3 &cameraPosition, float specularExponent, mat3 Camera_Orientation) {
    glm::vec3 interpolatedNormal = getInterpolatedNormal(intersection.intersectedTriangle, intersection.u, intersection.v);
    Colour colour = intersection.intersectedTriangle.colour;
    float totalProximityLighting = 0.0f;
    float totalIncidenceLighting = 0.0f;
    float totalSpecularLighting = 0.0f;

    for (const auto& lightPosition : lightSources) {
        float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
        float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, interpolatedNormal);
        float specularLighting = getSpecularIntensity(lightPosition, intersection.intersectionPoint, interpolatedNormal, specularExponent, cameraPosition, Camera_Orientation);

        totalProximityLighting += ProximityLighting;
        totalIncidenceLighting += IncidenceLighting;
        totalSpecularLighting += specularLighting;
    }

    // 平均化来自所有光源的照明
    float avgProximityLighting = totalProximityLighting / lightSources.size();
    float avgIncidenceLighting = totalIncidenceLighting / lightSources.size();
    float avgSpecularLighting = totalSpecularLighting / lightSources.size();

    colour.red = static_cast<uint8_t>(std::min(colour.red * avgProximityLighting * avgIncidenceLighting + avgSpecularLighting * 255, 255.0f));
    colour.green = static_cast<uint8_t>(std::min(colour.green * avgProximityLighting * avgIncidenceLighting + avgSpecularLighting * 255, 255.0f));
    colour.blue = static_cast<uint8_t>(std::min(colour.blue * avgProximityLighting * avgIncidenceLighting + avgSpecularLighting * 255, 255.0f));

    return colour;
}





glm::vec3 getRayDirection(float x,float y,glm::vec3 camera_position,mat3 Camera_Orientation){
    float scalingFactor = 150.0f;
    float focal_length = 2.0f;

    double x1 = (x- camera_position.x - (WIDTH / 2)) / (scalingFactor * focal_length);
    double y1 = -(y-camera_position.y - (HEIGHT / 2)) / (scalingFactor * focal_length);
    double z1 = -focal_length ;

    glm::vec3 ray_direction = glm::vec3(x1 , y1 , z1);
    glm::vec3 direction= normalize(Camera_Orientation*ray_direction);
    return direction;
}



Colour CalculateColor(const std::vector<ModelTriangle>& triangles, RayTriangleIntersection& intersection,  glm::vec3& camera_position,  vector<vec3>& lightSources, float specularExponent,  mat3& Camera_Orientation, bool diffuse, bool specular, bool gouraud, bool phong) {
    Colour colour;
    bool inShadow = IsShadow(intersection, lightSources, triangles);
    glm::vec3 normal = calculateNormal(intersection.intersectedTriangle, intersection);

    if (gouraud) {
        colour = gouraudShading(intersection, lightSources, camera_position, specularExponent, Camera_Orientation,triangles);
    } else if (phong) {
        colour = phongShading(intersection, lightSources, camera_position, specularExponent, Camera_Orientation);
    } else if (diffuse || specular) {
        float brightness = getBrightness(intersection, lightSources, normal, camera_position, specularExponent,triangles);
        //float SpecularLighting = getSpecularIntensity(lightSources, intersection.intersectionPoint, normal, specularExponent, camera_position, Camera_Orientation);

        if (diffuse) {
            colour = intersection.intersectedTriangle.colour;
            colour.red *= std::min(brightness, 1.0f);
            colour.green *= std::min(brightness, 1.0f);
            colour.blue *= std::min(brightness, 1.0f);
        }
        if (specular) {
            colour = intersection.intersectedTriangle.colour;
            brightness= getVertexBrightness(intersection,lightSources,normal,camera_position,256.0f,Camera_Orientation,triangles);
            colour.red *= std::min(brightness, 1.0f);
            colour.green *= std::min(brightness, 1.0f);
            colour.blue *= std::min(brightness, 1.0f);

//            colour.red = std::min(static_cast<int>(colour.red * brightness + SpecularLighting * 255), 255);
//            colour.green = std::min(static_cast<int>(colour.green * brightness + SpecularLighting * 255), 255);
//            colour.blue = std::min(static_cast<int>(colour.blue * brightness + SpecularLighting * 255), 255);
        }
    }

    if (inShadow) {
        colour.red /= 3;
        colour.green /= 3;
        colour.blue /= 3;
    }

    return colour;
}
void DrawRay(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation,vector<vec3>& lightSources,float specularExponent,
             bool diffuse,bool specular,bool gouraud,bool phong) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            glm::vec3 rayDirection= getRayDirection(x,y,camera_position,Camera_Orientation);
            RayTriangleIntersection intersection = Get_closest_intersection(camera_position, rayDirection, triangles);
            if (!isinf(intersection.distanceFromCamera)) {
                 Colour colour = CalculateColor(triangles,intersection, camera_position, lightSources, specularExponent, Camera_Orientation, diffuse, specular, gouraud, phong);

                    uint32_t color_value =
                            (255 << 24) + (static_cast<int>(colour.red) << 16) + (static_cast<int>(colour.green) << 8) +
                            static_cast<int>(colour.blue);
                    window.setPixelColour(x, y, color_value);
                }

            }
        }
    }
