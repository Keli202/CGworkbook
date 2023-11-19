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
#include <cmath>
//#include "Ray.h"
using namespace std;
using namespace glm;
#ifndef REDNOISE_CAMERA_H
#define REDNOISE_CAMERA_H
#endif //REDNOISE_CAMERA_H
mat3 Camera_Orientation(vec3( -1.0,    0.0,    0.0),vec3( 0.0, 1.0,0.0),vec3( 0.0, 0.0, 1.0));

vec3 center = vec3(0.0f, 0.0f, 0.0f);

mat3 rotation_y(float t) {
    return mat3(vec3( cos(t), 0.0, sin(t)),vec3(    0.0, 1.0,    0.0),vec3(-sin(t), 0.0, cos(t)));
}

mat3 rotation_x(float t) {
    return mat3(vec3( 1.0,    0.0,    0.0),vec3( 0.0, cos(t),-sin(t)),vec3( 0.0, sin(t), cos(t)));
}

vec3 lightPosition(0.0f, 0.6f, 0.8f);

mat3 lookAt(vec3 cameraPosition, vec3 lookPoint) {
    vec3 forward = normalize(cameraPosition-lookPoint);
    vec3 vertical = vec3(0.0f, 1.0f, 0.0f);
    vec3 right = normalize(cross(vertical, forward));
    vec3 up = normalize(cross(forward, right));
    return mat3(-right, up, forward);
}

void orbitAndLookAt(vec3& cameraPosition, vec3 orbitCenter, float angle) {
    vec3 relativePosition = cameraPosition-orbitCenter;

    float cosA = cos(angle);
    float sinA = sin(angle);

    vec3 newPosition;
    newPosition.x = relativePosition.x * cosA + relativePosition.z * sinA;
    newPosition.y = relativePosition.y;
    newPosition.z = -relativePosition.x * sinA + relativePosition.z * cosA;
    cameraPosition = orbitCenter + newPosition;

    Camera_Orientation = lookAt(cameraPosition, orbitCenter);
}

static bool orbitEnabled = false;
bool render=false;
bool ray=false;
bool wire=false;
bool light=false;

void changePosition(const std::vector<ModelTriangle>& modelTriangles,vec3& cameraPosition,SDL_Event event, DrawingWindow &window,std::vector<std::vector<float>> depthBuffer){
    float a=0.1f;
    float t=M_PI/180;
    float lightMoveStep = 0.1f;
    if (event.type == SDL_KEYDOWN) {

            if (event.key.keysym.sym == SDLK_a) { cameraPosition.x += a; }
            else if (event.key.keysym.sym == SDLK_d) { cameraPosition.x -= a; }
            else if (event.key.keysym.sym == SDLK_w) { cameraPosition.y -= a; }
            else if (event.key.keysym.sym == SDLK_s) { cameraPosition.y += a; }
            else if (event.key.keysym.sym == SDLK_o) { cameraPosition.z += a; }
            else if (event.key.keysym.sym == SDLK_p) { cameraPosition.z -= a; }
            else if (event.key.keysym.sym == SDLK_k) { cameraPosition = cameraPosition * rotation_x(-t); }
            else if (event.key.keysym.sym == SDLK_i) { cameraPosition = cameraPosition * rotation_x(t); }
            else if (event.key.keysym.sym == SDLK_l) { cameraPosition = cameraPosition * rotation_y(-t); }
            else if (event.key.keysym.sym == SDLK_j) { cameraPosition = cameraPosition * rotation_y(t); }
            else if (event.key.keysym.sym == SDLK_LEFT) { Camera_Orientation = Camera_Orientation * rotation_y(-t); }
            else if (event.key.keysym.sym == SDLK_RIGHT) { Camera_Orientation = Camera_Orientation * rotation_y(t); }
            else if (event.key.keysym.sym == SDLK_UP) { Camera_Orientation = Camera_Orientation * rotation_x(-t); }
            else if (event.key.keysym.sym == SDLK_DOWN) { Camera_Orientation = Camera_Orientation * rotation_x(t); }
            else if (event.key.keysym.sym == SDLK_z) { lightPosition.x -= lightMoveStep; }
            else if (event.key.keysym.sym == SDLK_x) { lightPosition.x += lightMoveStep; }
            else if (event.key.keysym.sym == SDLK_c) { lightPosition.y += lightMoveStep; }
            else if (event.key.keysym.sym == SDLK_v) { lightPosition.y -= lightMoveStep; }
            else if (event.key.keysym.sym == SDLK_b) { lightPosition.z -= lightMoveStep; }
            else if (event.key.keysym.sym == SDLK_n) { lightPosition.z += lightMoveStep; }
            else if (event.key.keysym.sym == SDLK_r) {
                float r = sqrt(cameraPosition.x * cameraPosition.x + cameraPosition.z * cameraPosition.z);
                float theta = atan2(cameraPosition.z, cameraPosition.x);
                float deltaTheta = 0.05f;
                theta += deltaTheta;
                cameraPosition.x = r * cos(theta);
                cameraPosition.z = r * sin(theta);
            } else if (event.key.keysym.sym == SDLK_g) {
                orbitEnabled = !orbitEnabled;
                if (orbitEnabled) {
                    orbitAndLookAt(cameraPosition, center, M_PI / 180);
                }
            } else if (event.key.keysym.sym == SDLK_1) {

                render = !render;
            } else if (event.key.keysym.sym == SDLK_2) {
                ray = !ray;
            } else if (event.key.keysym.sym == SDLK_3) {
                wire = !wire;
            }else if (event.key.keysym.sym == SDLK_4) {
                //light=!light;
                diffuse = !diffuse;
            } else if (event.key.keysym.sym == SDLK_5) {
                specular = !specular;
            } else if (event.key.keysym.sym == SDLK_6) {
                light = !light;
            } else if (event.key.keysym.sym == SDLK_7) {
                gouraud = !gouraud;
            } else if (event.key.keysym.sym == SDLK_8) {
                phong = !phong;
            }


    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }

    //cout<<lightPosition[0]<<","<<lightPosition[1]<<","<<lightPosition[2]<<endl;
    window.clearPixels();
    if(render){    RenderScene(window, modelTriangles,cameraPosition,Camera_Orientation, depthBuffer);}
    if(ray){draw_raytrace(modelTriangles, window, cameraPosition, Camera_Orientation);}
    if(wire) {RenderTriangle(window, modelTriangles,cameraPosition,Camera_Orientation, depthBuffer);}
    //if(light) { drawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition);}
    //drawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition);
    if(diffuse) { DrawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition,256.0f);}
    if(specular) { DrawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition,256.0f);}
    if(gouraud) { DrawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition,256.0f);}
    if(phong) { DrawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition,256.0f);}


    //if(light) { drawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition);}
    //if(diffuse) { drawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition);}
    //if(specular) { drawRaytrace(modelTriangles, window, cameraPosition, Camera_Orientation,lightPosition);}

}
