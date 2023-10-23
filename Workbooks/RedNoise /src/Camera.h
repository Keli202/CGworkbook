#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#define WIDTH 320
#define HEIGHT 240
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <ModelTriangle.h>

using namespace std;
using namespace glm;
#ifndef REDNOISE_CAMERA_H
#define REDNOISE_CAMERA_H
#endif //REDNOISE_CAMERA_H
void changePosition(const std::vector<ModelTriangle>& modelTriangles,vec3& cameraPosition,SDL_Event event, DrawingWindow &window){
    cout<<"c1:"<<cameraPosition.x<<endl;

    if (event.type == SDL_KEYDOWN) {
        cout<<"c3:"<<cameraPosition.x<<endl;
        if (event.key.keysym.sym == SDLK_a) { cameraPosition.x+=0.1f;}
        else if (event.key.keysym.sym == SDLK_d) {cameraPosition.x-=0.1f;}
        else if (event.key.keysym.sym == SDLK_w) {cameraPosition.y-=0.1f;}
        else if (event.key.keysym.sym == SDLK_s) {cameraPosition.y+=0.1f;}
    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    cout<<"c2:"<<cameraPosition.x<<endl;
    window.clearPixels();
    RenderScene(window, modelTriangles,cameraPosition);

}
//const std::vector<ModelTriangle> modelTriangles;
//vec3 cameraPosition;
//void handleEvent(SDL_Event event, DrawingWindow &window) {
//    if (event.type == SDL_KEYDOWN) {
//        if (event.key.keysym.sym == SDLK_LEFT) cameraPosition.x+=0.1;
//        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
//        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
//        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
//    }
//    else if (event.type == SDL_MOUSEBUTTONDOWN) {
//        window.savePPM("output.ppm");
//        window.saveBMP("output.bmp");
//    }
//    window.clearPixels();
//    RenderScene(window, modelTriangles);
//}