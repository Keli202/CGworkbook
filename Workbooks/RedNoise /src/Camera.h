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
using namespace std;
using namespace glm;
#ifndef REDNOISE_CAMERA_H
#define REDNOISE_CAMERA_H
#endif //REDNOISE_CAMERA_H
vec3 Camera_Up = vec3(0.0, 1.0, 0.0);
//mat3 Camera_Orientation(vec3( 1.0,    0.0,    0.0),vec3( 0.0, 1.0,0.0),vec3( 0.0, 0.0, 1.0));
mat3 Camera_Orientation(vec3( 1.0,    0.0,    0.0),Camera_Up,vec3( 0.0, 0.0, 1.0));

mat3 rotation_y(float t) {
    return mat3(vec3( cos(t), 0.0, sin(t)),vec3(    0.0, 1.0,    0.0),vec3(-sin(t), 0.0, cos(t)));
}
mat3 rotation_x(float t) {
    return mat3(vec3( 1.0,    0.0,    0.0),vec3( 0.0, cos(t),-sin(t)),vec3( 0.0, sin(t), cos(t)));
}

mat3 lookAt(vec3 cameraPosition, vec3 target) {
    vec3 forward = normalize(target - cameraPosition);
    vec3 vertical = vec3(0.0f, 1.0f, 0.0f);

    // 根据提示计算"Right"和"Up"向量
    vec3 right = normalize(cross(vertical, forward));
    vec3 up = normalize(cross(forward, right));

    return mat3(right, up, forward);
}


void orbitAndLookAt(vec3& cameraPosition, vec3 orbitCenter, float angleIncrement) {
    vec3 relativePos = cameraPosition - orbitCenter;

    float cosA = cos(angleIncrement);
    float sinA = sin(angleIncrement);

    vec3 newPosition;
    newPosition.x = relativePos.x * cosA - relativePos.z * sinA;
    newPosition.y = relativePos.y;
    newPosition.z = relativePos.x * sinA + relativePos.z * cosA;


    cameraPosition = orbitCenter + newPosition;


    mat3 orientation = lookAt(cameraPosition, orbitCenter);

}


void changePosition(const std::vector<ModelTriangle>& modelTriangles,vec3& cameraPosition,SDL_Event event, DrawingWindow &window,std::vector<std::vector<float>> depthBuffer){
    float a=0.1f;
    float t=M_PI/180;
    if (event.type == SDL_KEYDOWN) {

        if (event.key.keysym.sym == SDLK_a) { cameraPosition.x+=a;}
        else if (event.key.keysym.sym == SDLK_d) {cameraPosition.x-=a;}
        else if (event.key.keysym.sym == SDLK_w) {cameraPosition.y-=a;}
        else if (event.key.keysym.sym == SDLK_s) {cameraPosition.y+=a;}
        else if (event.key.keysym.sym == SDLK_o) {cameraPosition.z+=a;}
        else if (event.key.keysym.sym == SDLK_p) {cameraPosition.z-=a;}
        else if (event.key.keysym.sym == SDLK_k) {cameraPosition = cameraPosition*rotation_x(-t);}
        else if (event.key.keysym.sym == SDLK_i) {cameraPosition =cameraPosition*rotation_x(t);}
        else if (event.key.keysym.sym == SDLK_l) {cameraPosition = cameraPosition*rotation_y(-t);}
        else if (event.key.keysym.sym == SDLK_j) {cameraPosition = cameraPosition*rotation_y(t);}
        else if (event.key.keysym.sym == SDLK_LEFT) {Camera_Orientation=Camera_Orientation*rotation_y(-t);}
        else if (event.key.keysym.sym == SDLK_RIGHT) {Camera_Orientation=Camera_Orientation*rotation_y(t);}
        else if (event.key.keysym.sym == SDLK_UP) {Camera_Orientation=Camera_Orientation*rotation_x(-t);}
        else if (event.key.keysym.sym == SDLK_DOWN) {Camera_Orientation=Camera_Orientation*rotation_x(t);}
        else if (event.key.keysym.sym == SDLK_r) {float r = sqrt(cameraPosition.x*cameraPosition.x + cameraPosition.z*cameraPosition.z);
            float theta = atan2(cameraPosition.z, cameraPosition.x);
                float deltaTheta = 0.01f;
                theta += deltaTheta;

                cameraPosition.x = r * cos(theta);
                cameraPosition.z = r * sin(theta);

        }
//        else if (event.key.keysym.sym == SDLK_g){
//            vec3 center = vec3(0.0f, 0.0f, 0.0f);  // 场景的中心
//            float angleIncrement = M_PI / 180;  // 角度增量，例如每帧旋转1度
//            orbitAndLookAt(cameraPosition, center, angleIncrement);
//        }

    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    //cout<<"c2:"<<Camera_Orientation<<endl;
    window.clearPixels();
    RenderScene(window, modelTriangles,cameraPosition,Camera_Orientation, depthBuffer);

}
