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

std::vector<glm::vec3> lightSources() {
     glm::vec3 center(0.0f,0.6f,0.8f);
     float radius=0.1f;
     int numSources=1;
     float height=0.6f;
    std::vector<glm::vec3> lightSources;
    for (int i = 0; i < numSources; ++i) {
        float angle = 2.0f * M_PI * i / numSources;
        float x = center.x + radius * cos(angle);
        float z = center.z + radius * sin(angle);
        lightSources.push_back(glm::vec3(x, height, z));
    }
    return lightSources;
}



vec3 lightPosition(0.0f,0.6f,1.2f);
//vector <vec3> lightSource{
//
//        vec3(0.2f, 0.6f, 1.0f),
//        vec3(-0.2f, 0.6f, 1.0f),
//        vec3(0.2, 0.6, 1.2),
//         vec3(-0.2, 0.6, 1.2),
//         vec3(0.0, 0.6, 1.0),
//         vec3(0.0, 0.6, 1.2),
//         vec3(0.0, 0.6, 0.8),
//         vec3(-0.2, 0.6, 0.8),
//         vec3(0.1, 0.6, -0.1),
//         vec3(0.1, 0.6, 0.0),
////         vec3(0.1, 0.8, 0.1),
////         vec3(-0.1, 0.9, -0.1),
////         vec3(-0.1, 0.9, 0.0),
////         vec3(-0.1, 0.9, 0.1),
////         vec3(0.0, 0.9, -0.1),
////         vec3(0.0, 0.9, 0.0),
////         vec3(0.0, 0.9, 0.1),
////         vec3(0.1, 0.9, -0.1),
////         vec3(0.1, 0.9, 0.0),
////         vec3(0.1, 0.9, 0.1),
////        vec3(-0.2, 1.0, 1.8),
////        vec3(-0.2, 1.0, 2.0),
////        vec3(-0.2, 1.0, 2.2),
////        vec3(0.0, 1.0, 1.8),
////        vec3(0.0, 1.0, 2.0),
////        vec3(0.0, 1.0, 2.2),
////        vec3(0.2, 1.0, 1.8),
////        vec3(0.2, 1.0, 2.0),
////        vec3(0.2, 1.0, 2.2)
//
//};

TextureMap TextureFile("../texture_2.ppm");



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

vec3 refract(const vec3& incident, const vec3& normal, float indexOfRefraction) {
    float cosi = clamp(-1.0f, 1.0f, dot(incident, normal));
    float etai = 1, etat = indexOfRefraction;
    vec3 n = normal;
    if (cosi < 0) {
        cosi = -cosi;
    } else {
        std::swap(etai, etat);
        n= -normal;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(0,0,0) : eta * incident + (eta * cosi - sqrtf(k)) * n;
}

RayTriangleIntersection Get_closest_intersection(vec3 cameraPosition, vec3 rayDirection, const std::vector<ModelTriangle>& triangles,int depth) {
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
                closestIntersection.u = u;
                closestIntersection.v = v;
                closestIntersection=RayTriangleIntersection(intersectionPoint,closestIntersection.distanceFromCamera,triangle,closestIntersection.triangleIndex,u,v);

            }
        }
    }
    if(depth>=5){
        closestIntersection.intersectedTriangle.colour = Colour(255,255,255);
        return closestIntersection;
    }
    if(closestIntersection.intersectedTriangle.mirror && depth<5){

        vec3 normal= calculateFaceNormal(closestIntersection.intersectedTriangle);
        glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
        RayTriangleIntersection mirrorIntersection=Get_closest_intersection(closestIntersection.intersectionPoint + reflectedDirection * 0.001f,reflectedDirection,triangles,depth+1);
        if(!isinf(mirrorIntersection.distanceFromCamera)){
            closestIntersection=mirrorIntersection;
        }

    }else if(closestIntersection.intersectedTriangle.glass && depth<5){

        vec3 normal= calculateFaceNormal(closestIntersection.intersectedTriangle);
        float refractiveIndex = 1.5;

        glm::vec3 reflectedDirection = refract(rayDirection, normal, refractiveIndex);

        RayTriangleIntersection glassIntersection=Get_closest_intersection(closestIntersection.intersectionPoint + reflectedDirection * 0.001f,reflectedDirection,triangles,depth+1);
        if(!isinf(glassIntersection.distanceFromCamera)){
            closestIntersection=glassIntersection;
        }

    }
    else if(closestIntersection.intersectedTriangle.texture && depth<5){


    }

    return closestIntersection;
}


int IsShadow(RayTriangleIntersection intersection, vector<vec3>& lightSources, const vector<ModelTriangle>& triangles) {
    int unblockedSamples = 0;
    for (const auto& lightPosition : lightSources) {
        vec3 shadow_ray = normalize(lightPosition - intersection.intersectionPoint);
        bool blocked = false;
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
                    blocked = true;
                    break;
                }
            }
        }
        if (!blocked) {
            unblockedSamples++;
        }
    }
    return unblockedSamples;
}







float getProximityLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint) {
    float distance = glm::length(lightPosition - intersectionPoint);
    return std::min(1.0f, std::max(0.0f, 20.0f / static_cast<float>(4 * M_PI * distance * distance)));
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






float getBrightness( RayTriangleIntersection &intersection, glm::vec3 &normal, glm::vec3 cameraPosition, float specularExponent, const vector<ModelTriangle>& triangles, mat3 Camera_Orientation, bool diffuse, bool specular){
    float totalDiffuseIntensity = 0.0f;
    float maxSpecularIntensity = 0.0f;
    std::vector<glm::vec3> lightSource = lightSources();

    for (const auto& lightPosition : lightSource) {
        float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
        float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, normal);
        totalDiffuseIntensity += ProximityLighting * IncidenceLighting;
        if (specular) {
            float specularLighting = getSpecularIntensity(lightPosition, intersection.intersectionPoint, normal, specularExponent, cameraPosition, Camera_Orientation);
            maxSpecularIntensity = std::max(maxSpecularIntensity, specularLighting);
        }
    }

    float ambientIntensity = 0.1f;
    float totalBrightness;
    if(diffuse){
        totalBrightness = ambientIntensity + totalDiffuseIntensity / lightSource.size();
    }

    if (specular) {
        totalBrightness = ambientIntensity + totalDiffuseIntensity / lightSource.size()+maxSpecularIntensity;
    }

    return glm::clamp(totalBrightness, 0.2f, 1.0f);
}


float getVertexBrightness(const RayTriangleIntersection &intersection, vec3& lightPosition, const glm::vec3 &normal, const glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation,const vector<ModelTriangle>& triangles) {
    float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
    //incidence
    float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, normal);

    float specularLighting = getSpecularIntensity(lightPosition, intersection.intersectionPoint, normal,specularExponent, cameraPosition, Camera_Orientation);
    //get brightness
    float ambientIntensity = 0.1f;

    float brightness = ambientIntensity + ProximityLighting * IncidenceLighting + specularLighting;

    brightness = glm::clamp(brightness, 0.0f, 1.0f);

    return brightness;
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


Colour gouraudShading(const RayTriangleIntersection &intersection,  glm::vec3 &lightPosition, const glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation,const vector<ModelTriangle>& triangles) {
    const ModelTriangle &triangle = intersection.intersectedTriangle;
    std::vector<float> vertexLightIntensities(3);
    for (int i = 0; i < 3; ++i) {
        glm::vec3 vertex = triangle.vertices[i];
        glm::vec3 vertexNormal = triangle.vertexNormals[i];

        float vertexBrightness = getVertexBrightness(intersection, lightPosition, vertexNormal, cameraPosition, specularExponent,Camera_Orientation,triangles);
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

    return colour;
}


Colour phongShading(RayTriangleIntersection &intersection, glm::vec3 &lightPosition, glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation) {
    glm::vec3 interpolatedNormal = getInterpolatedNormal(intersection.intersectedTriangle, intersection.u, intersection.v);

    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);
    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,interpolatedNormal);
    float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,interpolatedNormal,specularExponent,cameraPosition,Camera_Orientation);

    Colour colour = intersection.intersectedTriangle.colour;
    colour.red = static_cast<uint8_t>(std::min(colour.red * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
    colour.green = static_cast<uint8_t>(std::min(colour.green * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
    colour.blue = static_cast<uint8_t>(std::min(colour.blue * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));

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




Colour CalculateColor(const vector<ModelTriangle>& triangles,RayTriangleIntersection &intersection, glm::vec3 &cameraPosition, float specularExponent, mat3 &Camera_Orientation,  glm::vec3 rayDirection) {
    Colour colour;
    vector<vec3> lightSource=lightSources();
    vector<ModelTriangle> sphere1;
    vector<ModelTriangle> sphere2;
    vector<ModelTriangle> sphere3;

    vector<ModelTriangle> normalTriangles;



    glm::vec3 normal = calculateNormal(intersection.intersectedTriangle, intersection);
    int unblockedSamples = IsShadow(intersection, lightSource, triangles);
    float shadowIntensity = static_cast<float>(unblockedSamples) / static_cast<float>(lightSource.size());

    for(auto& triangle:triangles){
        if(triangle.material=="sphere1"){
            sphere1.push_back(triangle);
        }else if(triangle.material=="sphere2"){
            sphere2.push_back(triangle);
        }else if(triangle.material=="sphere3"){
            sphere3.push_back(triangle);
        }else{
            normalTriangles.push_back(triangle);
        }
    }
    if(intersection.intersectedTriangle.material=="sphere3"){
        colour = gouraudShading(intersection, lightPosition, cameraPosition, specularExponent, Camera_Orientation, triangles);
        //cout<<"g:"<<intersection.intersectedTriangle.material<<endl;
        //cout<<"g:"<<colour<<endl;

    }
    else if(intersection.intersectedTriangle.material=="sphere1"){
        colour = phongShading(intersection, lightPosition, cameraPosition, specularExponent, Camera_Orientation);
        //cout<<"p:"<<intersection.intersectedTriangle.material<<endl;

    }
    else if(intersection.intersectedTriangle.material=="sphere1"){
        //specular
        float brightness = getBrightness(intersection, normal, cameraPosition, specularExponent,sphere3,Camera_Orientation,false,true);
            shadowIntensity=clamp(shadowIntensity,0.4f,1.0f);
            brightness *= shadowIntensity;
            colour = intersection.intersectedTriangle.colour;
            colour.red *= std::min(brightness, 1.0f);
            colour.green *= std::min(brightness, 1.0f);
            colour.blue *= std::min(brightness, 1.0f);
        //cout<<"g:"<<colour<<endl;

    }
    else  {
       float brightness = getBrightness(intersection, normal, cameraPosition, specularExponent, normalTriangles,Camera_Orientation, true, false);
       shadowIntensity = clamp(shadowIntensity, 0.4f, 1.0f);
       brightness *= shadowIntensity;
       colour = intersection.intersectedTriangle.colour;
       colour.red *= std::min(brightness, 1.0f);
       colour.green *= std::min(brightness, 1.0f);
       colour.blue *= std::min(brightness, 1.0f);
   }

        if (intersection.intersectedTriangle.mirror) {
            glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
            RayTriangleIntersection reflectedIntersection = Get_closest_intersection(intersection.intersectionPoint + reflectedDirection * 0.001f, reflectedDirection, triangles,0);
            if (!isinf(intersection.distanceFromCamera)) {
                Colour reflectedColour = CalculateColor(triangles, reflectedIntersection, cameraPosition,
                                                        specularExponent, Camera_Orientation, rayDirection);
                colour = reflectedColour;
            }

    }else if(intersection.intersectedTriangle.glass){

            float reflectivity = 0.1f;

            glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
            RayTriangleIntersection reflectedIntersection = Get_closest_intersection(intersection.intersectionPoint + reflectedDirection * 0.001f, reflectedDirection, triangles,0);

            vec3 refractedRay = refract(rayDirection, normal, 1.5);
            RayTriangleIntersection refractionIntersection = Get_closest_intersection(intersection.intersectionPoint + refractedRay * 0.001f, refractedRay, triangles, 0);
            Colour reflectedColour =reflectedIntersection.intersectedTriangle.colour;
            Colour refractionColour =refractionIntersection.intersectedTriangle.colour;
            if (!isinf(intersection.distanceFromCamera)) {
                 reflectedColour = CalculateColor(triangles, reflectedIntersection, cameraPosition,
                                                     specularExponent, Camera_Orientation, rayDirection);

                     refractionColour = CalculateColor(triangles, refractionIntersection, cameraPosition, specularExponent, Camera_Orientation, refractedRay);

            }
            colour.red = reflectedColour.red * reflectivity + refractionColour.red * (1 - reflectivity);
            colour.green = reflectedColour.green * reflectivity + refractionColour.green * (1 - reflectivity);
            colour.blue = reflectedColour.blue * reflectivity +refractionColour.blue * (1 - reflectivity);



        }
        //Texture
//        else if(intersection.intersectedTriangle.texture){
//            TextureMap& textureFile=TextureFile;

              //cout<<"1"<<endl;
//            float proportion0 = 1 - (intersection.u + intersection.v);
//            float proportion1 = intersection.u;
//            float proportion2 = intersection.v;
//            TexturePoint tp0 = intersection.intersectedTriangle.texturePoints[0];
//            TexturePoint tp1 = intersection.intersectedTriangle.texturePoints[1];
//            TexturePoint tp2 = intersection.intersectedTriangle.texturePoints[2];
//
//            TexturePoint texturePoint;
//            texturePoint.x = (proportion0 * tp0.x) + (proportion1 * tp1.x) + (proportion2 * tp2.x);
//            texturePoint.y = (proportion0 * tp0.y) + (proportion1 * tp1.y) + (proportion2 * tp2.y);


//            ModelTriangle t = intersection.intersectedTriangle;
//
//            float x = ((1 - intersection.u - intersection.v) * t.texturePoints[0].x + intersection.u * t.texturePoints[1].x + intersection.v * t.texturePoints[2].x);
//            float y = ((1 - intersection.u - intersection.v) * t.texturePoints[0].y + intersection.u * t.texturePoints[1].y +intersection.v * t.texturePoints[2].y);
//
//            x *= TextureFile.width;
//            y *= TextureFile.height;
//            uint32_t colour_val = TextureFile.pixels[round(x) + (TextureFile.width * round(y))];

           //uint32_t colour_val = TextureFile.pixels[round(texturePoint.x) + (TextureFile.width * round(texturePoint.y))];
//

//
//            const ModelTriangle& triangle = intersection.intersectedTriangle;
//            float w = 1.0f - intersection.u - intersection.v;
//
//            TexturePoint texturePoint;
//            texturePoint.x = w * triangle.texturePoints[0].x + intersection.u * triangle.texturePoints[1].x + intersection.v * triangle.texturePoints[2].x;
//            texturePoint.y = w * triangle.texturePoints[0].y + intersection.u * triangle.texturePoints[1].y + intersection.v * triangle.texturePoints[2].y;
//
//            int x = static_cast<int>(texturePoint.x * TextureFile.width) %TextureFile.width;
//            int y = static_cast<int>(texturePoint.y * TextureFile.height) % TextureFile.height;
//            uint32_t colour_val =  TextureFile.pixels[y * TextureFile.width + x];

//
//            int red = (colour_val >> 16) & 0xFF;
//            int green = (colour_val >> 8) & 0xFF;
//            int blue = colour_val & 0xFF;
//            colour.red = red;
//            colour.green = green;
//            colour.blue = blue;






//            float u = intersection.u;
//            float v = intersection.v;
//            float w = 1.0f - u - v; // The third barycentric coordinate
//            TexturePoint tp0 = intersection.intersectedTriangle.texturePoints[0];
//            TexturePoint tp1 = intersection.intersectedTriangle.texturePoints[1];
//            TexturePoint tp2 = intersection.intersectedTriangle.texturePoints[2];
//            float textureX = (w * tp0.x + u * tp1.x + v * tp2.x) * (textureFile.width - 1);
//            float textureY = (w * tp0.y + u * tp1.y + v * tp2.y) * (textureFile.height - 1);
//
//            // Ensure we don't go out of bounds
//            textureX = fmod(textureX, textureFile.width);
//            textureY = fmod(textureY, textureFile.height);
//
//            // Sample the texture color
//            uint32_t colour_val = textureFile.pixels[static_cast<int>(textureY) * textureFile.width + static_cast<int>(textureX)];
//            int red = (colour_val >> 16) & 0xFF;
//            int green = (colour_val >> 8) & 0xFF;
//            int blue = colour_val & 0xFF;
//            colour.red = red;
//            colour.green = green;
//            colour.blue = blue;
//
//
//
//        }


    return colour;
}

void DrawRay(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation,float specularExponent) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            glm::vec3 rayDirection= getRayDirection(x,y,camera_position,Camera_Orientation);
            RayTriangleIntersection intersection = Get_closest_intersection(camera_position, rayDirection, triangles,1);
            if (!isinf(intersection.distanceFromCamera)) {
                 Colour colour = CalculateColor(triangles,intersection, camera_position, specularExponent, Camera_Orientation, rayDirection);

                    uint32_t color_value =
                            (255 << 24) + (static_cast<int>(colour.red) << 16) + (static_cast<int>(colour.green) << 8) +
                            static_cast<int>(colour.blue);
                    window.setPixelColour(x, y, color_value);
                }

            }
        }
    }
