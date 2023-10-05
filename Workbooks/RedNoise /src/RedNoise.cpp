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
using namespace std;
using namespace glm;
std::vector<float> interpolateSingleFloats(float from,float to, int numberOfValues){
    std::vector<float> results;
    float step = (to-from)/(numberOfValues-1);
    /* written by computer*/
    results.reserve(numberOfValues);
    for (int i=0; i<numberOfValues;i++){
        results.push_back(from+step*i);
    }

    return results;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from,glm::vec3 to, int numberOfValues){
    std::vector<glm::vec3> results;
    float stepX =(to.x-from.x)/(numberOfValues-1);
    float stepY =(to.y-from.y)/(numberOfValues-1);
    float stepZ =(to.z-from.z)/(numberOfValues-1);
    for (int i=0; i<numberOfValues;i++){
        glm::vec3 everyResult(from.x+stepX*i,from.y+stepY*i,from.z+stepZ*i);
        results.push_back(everyResult);
    }

    return results;

}
void drawColor(DrawingWindow &window){
    window.clearPixels();
    int Width = static_cast<int>(window.width);
    int Height = static_cast<int>(window.height);
    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    vector<vec3> Left = interpolateThreeElementValues(topLeft,bottomLeft,Height);
    vector<vec3> Right =interpolateThreeElementValues(topRight,bottomRight,Height);
    for (size_t y = 0; y < window.height; y++) {
        for (size_t x = 0; x < window.width; x++) {
            vec3 LeftValue = Left[y];
            vec3 RightValue = Right[y];
            vector<vec3> Row = interpolateThreeElementValues(LeftValue,RightValue,Width);
            vec3 colorValue = Row[x];
            uint32_t colour = (255<<24)+ (int(colorValue.x) << 16) + (int(colorValue.y) << 8) + int(colorValue.z);
            window.setPixelColour(x, y, colour);
        }
    }

}

void draw(DrawingWindow &window) {
    window.clearPixels();
    int Width = static_cast<int>(window.width);
    std::vector<float> gradient = interpolateSingleFloats(255.0f, 0.0f, Width);

    for (size_t y = 0; y < window.height; y++) {
        for (size_t x = 0; x < window.width; x++) {
            //float red = rand() % 256;
            //float green = 0.0;
            //float blue = 0.0;
            float greyValue = gradient[x];
            //uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            uint32_t colour = (255<<24)+ (int(greyValue) << 16) + (int(greyValue) << 8) + int(greyValue);
            window.setPixelColour(x, y, colour);
        }
    }
}

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



void drawLine(CanvasPoint from, CanvasPoint to, Colour inputColour, DrawingWindow &window){
    float xdistance = to.x - from.x;
    float ydistance = to.y - from.y;
    //abs() will calculate the absolute value of a number
    float numberOfSteps = std::max(abs(xdistance),abs(ydistance));
    //According to the distance to long, and then long distance one by one, you can see ppt.
    float xStepSize = xdistance/numberOfSteps;
    float yStepSize = ydistance/numberOfSteps;
    for (float i= 0.0; i<numberOfSteps;i++){
        float x =from.x+(i*xStepSize);
        float y =from.y+(i*yStepSize);
        uint32_t Colour = (255 << 24) + (inputColour.red << 16) + (inputColour.green << 8) + inputColour.blue;
        //round(x) is a mathematical function used to round a floating-point number x to the nearest integer
        window.setPixelColour(round(x), round(y), Colour);
    }
}


int main(int argc, char *argv[]) {
    std::vector<float> result;
    result = interpolateSingleFloats(2.2, 8.5, 7);
    for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
    std::cout << std::endl;

    glm::vec3 from(1.0, 4.0, 9.2);
    glm::vec3 to(4.0, 1.0, 9.8);
    std::vector<glm::vec3> ThreedResults= interpolateThreeElementValues(from,to,4);
    for (glm::vec3 vec : ThreedResults) {
        std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
    }
    //for(size_t i=0; i<ThreedResults.size(); i++) std::cout << ThreedResults[i] << " ";
    //std::cout << std::endl;

    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;
    while (true) {
        CanvasPoint from= CanvasPoint(0,0);
        CanvasPoint to = CanvasPoint(window.width/2,window.height/2);
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window);
        //draw(window);
        //drawColor(window);
        drawLine(from,to,Colour{25,234,34},window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }


}