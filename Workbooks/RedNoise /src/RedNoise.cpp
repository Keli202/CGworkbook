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

#include "drawRender.h"
#include "drawTriangle.h"
#include "light.h"
#include "testLight.h"
#include "Ray.h"
#include "Camera.h"
#include "draw.h"
#include "readOBJ.h"

using namespace std;
using namespace glm;


void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}



//week4
std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, std::numeric_limits<float>::infinity()));

void ClearDepthBuffer() {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            depthBuffer[x][y] = std::numeric_limits<float>::infinity();
        }
    }
}


void printMat3(const mat3& mat) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
}


int main(int argc, char *argv[]) {
   // week1();
    //loadOBJ("cornell-box.obj", 0.35);
//    for (const glm::vec3& vertex : modelVertices) {
//        std::cout << "Vertex: (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
//    }

    //std::cout << "Canvas Point: (" << canvasPoint.x << ", " << canvasPoint.y << ")" << std::endl;
//    for (int x = 0; x < WIDTH; x++) {
//        for (int y = 0; y < HEIGHT; y++) {
//            depthBuffer[x][y] = 0.0f;
//        }
//    }

    glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
     //std::vector<ModelTriangle> modelTriangles=readOBJ("../cornell-box.obj", 0.35);
     //std::vector<ModelTriangle> modelTriangles=readOBJ("../textured-cornell-box.obj", 0.35);
     std::vector<ModelTriangle> modelTriangles=readOBJ("../sphere.obj", 0.35);


    float focal_length = 2.0f;

    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    //mat3 Camera_Orientation(vec3( 1.0,    0.0,    0.0),vec3( 0.0, 1.0,0.0),vec3( 0.0, 0.0, 1.0));
    SDL_Event event;
    while (true) {
        window.clearPixels();
        ClearDepthBuffer();
        // We MUST poll for events - otherwise the window will freeze !
        //if (window.pollForInputEvents(event)) handleEvent(event, window);
        //week3(event,window);
       //week4(modelTriangles,window);
        //RenderScene(window, modelTriangles,cameraPosition);

        //draw_raytrace(modelTriangles, window, cameraPosition);
        modelTriangles= getVertexNormals(modelTriangles);
        changePosition(modelTriangles,cameraPosition,event,window,depthBuffer);

        //RenderScene(window, modelTriangles, cameraPosition, Camera_Orientation,depthBuffer);
        if(window.pollForInputEvents(event)) {
            handleEvent(event, window);
            //changePosition(modelTriangles,cameraPosition,event,window,depthBuffer,Camera_Orientation);

        }
        if (orbitEnabled) {
            float angle = -M_PI / 180;
            orbitAndLookAt(cameraPosition, center, angle);
           // cout<<"camPO:"<<cameraPosition[0]<<","<<cameraPosition[1]<<","<<cameraPosition[2]<<endl;
        }
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }


}