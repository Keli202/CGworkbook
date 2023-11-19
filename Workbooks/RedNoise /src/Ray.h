//
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
#ifndef REDNOISE_RAY_H
#define REDNOISE_RAY_H

#endif //REDNOISE_RAY_H

RayTriangleIntersection get_closest_intersection(vec3 cameraPosition, vec3 rayDirection, const std::vector<ModelTriangle>& triangles) {
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
                closestIntersection=RayTriangleIntersection(intersectionPoint,closestIntersection.distanceFromCamera,triangle,closestIntersection.triangleIndex,u,v);

            }
        }
    }

    return closestIntersection;
}
bool isShadow(RayTriangleIntersection intersection, vec3 lightPosition,const vector<ModelTriangle>& triangles){
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



    void draw_raytrace(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            float scalingFactor = 150.0f;
            float focal_length = 2.0f;

//            double x1 = -(x - (WIDTH / 2)) / (scalingFactor * focal_length)- camera_position.x;
//            double y1 = -(y - (HEIGHT / 2)) / (scalingFactor * focal_length)-camera_position.y;
//            double z1 = -focal_length - camera_position.z;
            double x1 = -(x- camera_position.x - (WIDTH / 2)) / (scalingFactor * focal_length);
            double y1 = -(y-camera_position.y - (HEIGHT / 2)) / (scalingFactor * focal_length);
            double z1 = -focal_length ;

            glm::vec3 ray_direction = glm::vec3(x1 , y1 , z1);
            glm::vec3 direction= normalize(Camera_Orientation*ray_direction);



            // Find the closest intersection
            RayTriangleIntersection intersection = get_closest_intersection(camera_position, direction, triangles);
            //cout<<"1:"<<rt_int<<endl;


            vec3 lightPosition(0.0f, 0.5f, 0.5f);
            if (!isinf(intersection.distanceFromCamera)) {
                bool inShadow = isShadow(intersection, lightPosition, triangles);
                Colour colour = intersection.intersectedTriangle.colour;
                if (inShadow) {
                    //cout<<"1"<<endl;
                    //cout<<"1"<<colour<<endl;
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
