//#include <CanvasTriangle.h>
//#include <DrawingWindow.h>
//#include <Utils.h>
//#include <fstream>
//#include <vector>
//#include <glm/glm.hpp>
//#define WIDTH 320*2.5
//#define HEIGHT 240*2.5
//#include <CanvasPoint.h>
//#include <Colour.h>
//#include <CanvasTriangle.h>
//#include <TextureMap.h>
//#include <ModelTriangle.h>
//
//#include "drawRender.h"
//#include "drawTriangle.h"
//#include "testLight.h"
//#include "light.h"
//#include "Ray.h"
//
//#include "Camera.h"
//#include "draw.h"
//#include "readOBJ.h"
//
//using namespace std;
//using namespace glm;
































//std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, std::numeric_limits<float>::infinity()));
//
//void ClearDepthBuffer() {
//    for (int x = 0; x < WIDTH; x++) {
//        for (int y = 0; y < HEIGHT; y++) {
//            depthBuffer[x][y] = std::numeric_limits<float>::infinity();
//        }
//    }
//}
//
//int main(int argc, char *argv[]) {
//
//
//    glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
//    //std::vector<ModelTriangle> modelTriangles=readOBJ("../cornell-box.obj", 0.35);
//    //std::vector<ModelTriangle> modelTriangles=readOBJ("../textured-cornell-box.obj", 0.35);
//    std::vector<ModelTriangle> modelTriangles=readOBJ("../sphere.obj", 0.35);
//
//
//
//    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
//    SDL_Event event;
//    while (true) {
//        window.clearPixels();
//        ClearDepthBuffer();
//
//        changePosition(modelTriangles,cameraPosition,event,window,depthBuffer);
//
//        if(window.pollForInputEvents(event)) {
//
//        }
//        if (orbitEnabled) {
//            float angle = -M_PI / 180;
//            orbitAndLookAt(cameraPosition, center, angle);
//        }
//        window.renderFrame();
//    }
//
//
//}