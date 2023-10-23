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
#include "drawRender.h"
#include "drawTriangle.h"
#include "Camera.h"
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



//week4
//std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, -numeric_limits<float>::infinity()));
//std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0.0f));

//void ClearDepthBuffer() {
//    for (int x = 0; x < WIDTH; x++) {
//        for (int y = 0; y < HEIGHT; y++) {
//            depthBuffer[x][y] = -std::numeric_limits<float>::infinity();
//        }
//    }
//}

void readOBJColour(const std::string& filename,std::vector<std::pair<std::string, Colour>>& palette) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open material file " << filename << std::endl;
        return;
    }
    string line;
    //vector<vec3> vertices;
    string content;
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0]=="newmtl"){
            //content +=line+ '\n';
            //palette.push_back(std::make_pair(currentMaterial, currentColour));
            content = tokens[1];
        }
        if(tokens[0]=="Kd"){
            float r, g, b;
            r = std::stof(tokens[1]) ;
            g = std::stof(tokens[2]) ;
            b = std::stof(tokens[3]) ;
            int red = static_cast<int>(r * 255);
            int green = static_cast<int>(g * 255);
            int blue = static_cast<int>(b * 255);
            palette.emplace_back(content, Colour(red, green, blue));
        }

    }
    file.close();
}

vector<ModelTriangle> readOBJ(const std::string& filename, float scale){
       ifstream file(filename);
       vector<ModelTriangle> triangles;
      vector<vec3> vertices;

       if(!file.is_open()){
           cerr<<"Error: Could not open file "<<filename<<std::endl;
           return triangles;
       }

        string line;

   // std::vector<Colour> colours;
    vector<pair<string, Colour>> palette;

    // Read the MTL file and populate the palette
    readOBJColour("../cornell-box.mtl",palette);
    string currentMaterial;

        while(getline(file,line)){
           if (line.empty()) {
               continue;
           }
           std::vector<std::string> tokens = split(line, ' ');
           if (tokens[0] == "v") {
               vec3 vertex;
               //cout<<tokens.size()<<endl;
               vertex.x= std::stof(tokens[1]);
               vertex.y = std::stof(tokens[2]);
               vertex.z = std::stof(tokens[3]);
               vertex *= scale;
               vertices.push_back(vertex);
           }
           if (tokens[0] == "f") {
               std::vector<std::string> parts = split(line, '/');
               int v1, v2, v3;
               v1 = std::stoi(tokens[1]) - 1;
               v2 = std::stoi(tokens[2]) - 1;
               v3 = std::stoi(tokens[3]) - 1;

               Colour triangleColour; // Default colour
               for (const auto& data : palette) {
                   if (data.first == currentMaterial) {
                       triangleColour = data.second;
                       //cout<<currentMaterial<<endl;
                       break;
                   }
               }
               //cout<<triangleColour<<endl;
               ModelTriangle triangle(vertices[v1], vertices[v2], vertices[v3], triangleColour);
               triangles.push_back(triangle);
//               ModelTriangle triangle(vertices[v1], vertices[v2], vertices[v3],Colour{255,255,255});
//               triangles.push_back(triangle);
           }
           if(tokens[0] == "usemtl"){
               currentMaterial = tokens[1];
           }
       }

       for (const ModelTriangle& triangle : triangles) {
            std::cout << triangle << std::endl;
        }

        file.close();
return triangles;
//return vertices;
}


//void UpdateDepthBuffer(const CanvasPoint& point) {
//    int x = static_cast<int>(point.x);
//    int y = static_cast<int>(point.y);
//   // glm::vec3 relativePosition = vertexPosition - cameraPosition;
//   // float depth = 1.0f / relativePosition.z;
//    if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT ) {
//        depthBuffer[x][y] = std::max(depthBuffer[x][y], 1.0f/point.depth);
//    }
//}
//void ClearDepthBuffer() {
//    for (int x = 0; x < WIDTH; x++) {
//        for (int y = 0; y < HEIGHT; y++) {
//            depthBuffer[x][y] = numeric_limits<float>::infinity();
//        }
//    }
//}




void week1(){
    //week1
    std::vector<float> result;
    result = interpolateSingleFloats(2.2, 8.5, 7);
    for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
    std::cout << std::endl;

 //week1
    glm::vec3 from(1.0, 4.0, 9.2);
    glm::vec3 to(4.0, 1.0, 9.8);
    std::vector<glm::vec3> ThreedResults= interpolateThreeElementValues(from,to,4);
    for (glm::vec3 vec : ThreedResults) {
        std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
    }
}
void week3(SDL_Event event, DrawingWindow &window){
    //draw(window);
    //drawColor(window);

    // //drawLine
    //drawLine(CanvasPoint(0,0),CanvasPoint(window.width/2,window.height/2),Colour{255,255,255},window);
    //drawLine(CanvasPoint(window.width-1,0),CanvasPoint(window.width/2,window.height/2),Colour{255,255,255},window);
    //drawLine(CanvasPoint(window.width / 2, 0),CanvasPoint(window.width / 2, window.height),Colour{255,255,255},window);
    //drawLine(CanvasPoint(window.width / 3, window.height / 2),CanvasPoint(2 * window.width / 3, window.height / 2),Colour{255,255,255},window);

    Colour randomColour(rand() % 256, rand() % 256, rand() % 256);
    drawKeyTriangle(event,randomColour,window);
//
    //filled triangle
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_f) {
            // Generate random triangle vertices and colors
            CanvasPoint vertex1(rand() % window.width, rand() % window.height);
            CanvasPoint vertex2(rand() % window.width, rand() % window.height);
            CanvasPoint vertex3(rand() % window.width, rand() % window.height);

            // Create the CanvasTriangle object and draw the filled triangle
            CanvasTriangle triangle({vertex1, vertex2, vertex3});
            drawFilledTriangles(triangle, Colour(rand() % 256, rand() % 256, rand() % 256),
                                Colour(255, 255, 255), window);
        } else if (event.type == SDL_MOUSEBUTTONDOWN) {
            window.savePPM("output.ppm");
            window.saveBMP("output.bmp");
        }
    }

        //fill textureMap
        TextureMap texture=TextureMap("../texture.ppm");
        // Define the vertices of the triangle
        CanvasPoint p1(160, 10);
        CanvasPoint p2(300, 230);
        CanvasPoint p3(10, 150);
        // Create a textured triangle
        p1.texturePoint =TexturePoint(195,5);
        p2.texturePoint =TexturePoint(395,380);
        p3.texturePoint =TexturePoint(65,330);
        CanvasTriangle triangle(p1, p2, p3);
        // Draw the textured triangle
        drawTexturedTriangle(triangle, texture, window);

}
//void week4(std::vector<ModelTriangle> modelTriangles,DrawingWindow &window){
//    window.clearPixels();
//   // ClearDepthBuffer();
//    glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
//    float focalLength = 2.0f;
//    // Image plane scaling factor
//    float scalingFactor = 150.0f;
//    // Load the vertex data for the Cornell Box model
//
//    CanvasTriangle canvasTriangle;
//    for(auto & modelTriangle : modelTriangles){
//        for (int i=0;i<3;i++) {
//            vec3 vertexPosition =modelTriangle.vertices[i];
//            glm::vec3 canvasPoint = getCanvasIntersectionPoint(cameraPosition, vertexPosition, focalLength,scalingFactor);
//            // draw white pixels
////                uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + 255;
////                window.setPixelColour(static_cast<int>(canvasPoint.x), static_cast<int>(canvasPoint.y), colour);
//            if (canvasPoint.x != -1 && canvasPoint.y != -1) {
//                CanvasPoint p = {canvasPoint.x, canvasPoint.y,canvasPoint.z};
//                //canvasTriangle.vertices[i].depth=1.0f/(vertexPosition.z-cameraPosition.z);
//                canvasTriangle.vertices[i] = p;
//                cout<<"p:"<<p<<endl;
//            }
//
//
//         //   UpdateDepthBuffer(canvasTriangle.vertices[i]);
//
////            float avgDepth = (canvasTriangle.vertices[i].depth / 1.0f);
////            cout<<"avgDepth1:"<<avgDepth<<endl;
////            cout<<"depth:"<<canvasTriangle.vertices[i].depth<<endl;
//        }
//         //   drawSpecialTriangle(canvasTriangle, modelTriangle.colour, window);
//           // drawFilledTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window);
//        float avgDepth = (canvasTriangle.vertices[0].depth + canvasTriangle.vertices[1].depth + canvasTriangle.vertices[2].depth) / 3.0f;
//        cout<<"avgDepth2:"<<avgDepth<<endl;
//
//        for(int i=0;i<3;i++){
////        if (avgDepth < depthBuffer[static_cast<int>(canvasTriangle.vertices[i].x)][static_cast<int>(canvasTriangle.vertices[i].y)]) {
////            drawSpecialTriangle(canvasTriangle, modelTriangle.colour, window);
////            //drawFilledTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window);
////        }
//
//        }
//
//    }
//}




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
    std::vector<ModelTriangle> modelTriangles=readOBJ("../cornell-box.obj", 0.35);
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;
    while (true) {
       // ClearDepthBuffer();
        // We MUST poll for events - otherwise the window will freeze !
        //if (window.pollForInputEvents(event)) handleEvent(event, window);
        //week3(event,window);


       //week4(modelTriangles,window);
        //RenderScene(window, modelTriangles,cameraPosition);

        changePosition(modelTriangles,cameraPosition,event,window);
        if(window.pollForInputEvents(event)) handleEvent(event,window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }


}