
//

#ifndef REDNOISE_DRAWRENDER_H
#define REDNOISE_DRAWRENDER_H

#endif //REDNOISE_DRAWRENDER_H
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

using namespace std;
using namespace glm;
//std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, (std::numeric_limits<float>::infinity())));
bool isTriangleInsideViewFrustrum(const CanvasTriangle &triangle) {
    for(int i = 0; i < 3; i++) {
        if(triangle.vertices[i].x < 0 || triangle.vertices[i].x > WIDTH || triangle.vertices[i].y < 0 || triangle.vertices[i].y > HEIGHT) {
            return false;
        }
    }
    return true;
}
vec3 getCanvasIntersectionPoint(vec3 cameraPosition,vec3 vertexPosition,float focalLength, float scalingFactor,mat3 Camera_Orientation){

    glm::vec3 relativePosition = vertexPosition - cameraPosition;
    //relativePosition =  Camera_Orientation*relativePosition+cameraPosition;
    relativePosition =  Camera_Orientation*relativePosition;
    // Calculate the projection coordinates of the vertices on the image plane
    double u = (-focalLength * relativePosition.x / relativePosition.z)*scalingFactor + (WIDTH / 2);
    double v = (focalLength * relativePosition.y / relativePosition.z)*scalingFactor + (HEIGHT / 2);
    if (u >= 0 && u < WIDTH && v >= 0 && v < HEIGHT) {
        float depth = 1.0f/(relativePosition.z);
        return vec3(u, v,depth);
    }
    else {
        return vec3(-1, -1,-1);
        //return vec3(-u, -v,-relativePosition.z);
        cout<<"vec:"<<-u<<endl;
    }

}

void drawRenderLine(CanvasPoint from, CanvasPoint to, Colour inputColour, DrawingWindow &window,vector<vector<float>>& depthBuffer){
    float xdistance = to.x - from.x;
    float ydistance = to.y - from.y;
    float zdistance = to.depth-from.depth;
    //abs() will calculate the absolute value of a number
    //float numberOfSteps = std::max(std::max(abs(xdistance),abs(ydistance)),abs(zdistance));
   float numberOfSteps =std::max(abs(xdistance),abs(ydistance));

    float xStepSize = xdistance/numberOfSteps;
    float yStepSize = ydistance/numberOfSteps;
    float zStepSize = zdistance/numberOfSteps;
    for (float i= 0; i<numberOfSteps;i++){
        float x =from.x+(i*xStepSize);
        float y =from.y+(i*yStepSize);
        float z =from.depth+(i*zStepSize);
        //cout<<"z:"<<from.depth<<endl;
        //float z=from.depth;
        cout<<"db0:"<<depthBuffer[round(x)][round(y)]<<endl;
        //cout<<"db0JJJJ:"<<depthBuffer[x][y]<<endl;
        cout<<"z0:"<<z<<endl;
        //if (fabs(z) < 1e-6) continue;
        //if (-1.0f/z>depthBuffer[round(x)][round(y)]){
            if (z<depthBuffer[round(x)][round(y)]){
           // cout<<"dbHHHHH:"<<depthBuffer[round(x)][round(y)]<<endl;
            //cout<<"zHHHH:"<<z<<endl;
            uint32_t Colour = (255 << 24) + (inputColour.red << 16) + (inputColour.green << 8) + inputColour.blue;
            //round(x) is a mathematical function used to round a floating-point number x to the nearest integer
            window.setPixelColour(round(x), round(y), Colour);
            depthBuffer[round(x)][round(y)]=z;
            // cout<<"db2:"<<depthBuffer[round(x)][round(y)]<<endl;
            // cout<<"z2:"<<z<<endl;
        }
    }
}


//void drawRenderTriangles(CanvasTriangle triangle, Colour fillColour, Colour LineColour, DrawingWindow &window,vector<vector<float>>& depthBuffer) {
//    CanvasPoint top = triangle[0];
//    CanvasPoint middle = triangle[1];
//    CanvasPoint bottom = triangle[2];
//    // Sort vertices by y-coordinate
//    if (top.y > middle.y) swap(top, middle);
//    if (top.y > bottom.y) swap(top, bottom);
//    if (middle.y > bottom.y) swap(middle, bottom);
//    // Calculate slopes for the top and bottom edges
//    float slope1 = (middle.x - top.x) / (middle.y - top.y);
//    float slope2 = (bottom.x - top.x) / (bottom.y - top.y);
//
//    float depthSlope1 = (middle.depth - top.depth) / (middle.y - top.y);
//    float depthSlope2 = (bottom.depth - top.depth) / (bottom.y - top.y);
//
//
//
//    float leftX = top.x;
//    float rightX = top.x;
//    float leftDepth = top.depth;
//    float rightDepth = top.depth;
//
//    for (float y=top.y;y<middle.y;y++){
//        CanvasPoint from(round(leftX),round(y),round(leftDepth));
//        CanvasPoint to(round(rightX),round(y),round(rightDepth));
//        drawRenderLine(from,to,fillColour,window,depthBuffer);
//        leftX+=slope1;
//        rightX+=slope2;
//        leftDepth+=depthSlope1;
//        rightDepth+=depthSlope2;
//
//    }
//     float leftX2=middle.x;
//     float leftDepth2=middle.depth;
//    float slope3 = (bottom.x - middle.x) / (bottom.y - middle.y);
//    float depthSlope3 = (bottom.depth - middle.depth)/(bottom.y - middle.y);
//    for (float y = middle.y ; y <=bottom.y; y++) {
//        CanvasPoint from(round(leftX2),round(y),round(leftDepth2)) ;
//        CanvasPoint to (round(rightX),round(y),round(rightDepth)) ;
//        drawRenderLine(from,to,fillColour,window,depthBuffer);
//        leftX2 += slope3;
//        rightX += slope2;
//        leftDepth+=depthSlope3;
//        rightDepth+=depthSlope2;
//    }
//
//    drawRenderLine(top, middle, LineColour, window,depthBuffer);
//    drawRenderLine(middle, bottom, LineColour, window,depthBuffer);
//    drawRenderLine(bottom, top, LineColour, window,depthBuffer);
//
//}
void drawRenderTriangles(CanvasTriangle triangle, Colour fillColour, Colour LineColour, DrawingWindow &window, vector<vector<float>>& depthBuffer) {

    if (triangle[0].y > triangle[1].y) std::swap(triangle[0], triangle[1]);
    if (triangle[0].y > triangle[2].y) std::swap(triangle[0], triangle[2]);
    if (triangle[1].y > triangle[2].y) std::swap(triangle[1], triangle[2]);

    float invSlope1 = (triangle[1].x - triangle[0].x) / (triangle[1].y - triangle[0].y);
    float invSlope2 = (triangle[2].x - triangle[0].x) / (triangle[2].y - triangle[0].y);
    float depthSlope1 = (triangle[1].depth - triangle[0].depth) / (triangle[1].y - triangle[0].y);
    float depthSlope2 = (triangle[2].depth - triangle[0].depth) / (triangle[2].y - triangle[0].y);

    float x1 = triangle[0].x;
    float x2 = triangle[0].x;
    float depth1 = triangle[0].depth;
    float depth2 = triangle[0].depth;
        for (float y = triangle[0].y; y < triangle[1].y; y++) {
            float startX = std::min(x1, x2);
            float endX = std::max(x1, x2);
            float startDepth = std::min(depth1, depth2);

            for (float x = round(startX); x <= round(endX); x++) {
                if (x >= 0 && x < window.width && y >= 0 && y < window.height) {

                    //if (-1.0f / startDepth > depthBuffer[round(x)][round(y)]) {
                        if (startDepth <depthBuffer[round(x)][round(y)]) {
                        depthBuffer[round(x)][round(y)] = startDepth;
                        uint32_t Colour =
                                (255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
                       // if (isTriangleInsideViewFrustrum(triangle)) {
                        window.setPixelColour(round(x), round(y), Colour);
                        //  }else{break;}
                    }
                }
            }
            x1 += invSlope1;
            x2 += invSlope2;
            depth1 += depthSlope1;
            depth2 += depthSlope2;
        }

        float invSlope3 = (triangle[2].x - triangle[1].x) / (triangle[2].y - triangle[1].y);
        float depthSlope3 = (triangle[2].depth - triangle[1].depth) / (triangle[2].y - triangle[1].y);

        x1 = triangle[1].x;
        depth1 = triangle[1].depth;

        for (float y = triangle[1].y; y <= triangle[2].y; y++) {
            float startX = std::min(x1, x2);
            float endX = std::max(x1, x2);
            float startDepth = std::min(depth1, depth2);

            for (float x = round(startX); x <= round(endX); x++) {
                if (x >= 0 && x < window.width && y >= 0 && y < window.height) {

                    //if (-1.0f / startDepth > depthBuffer[round(x)][round(y)]) {
                        if ( startDepth < depthBuffer[round(x)][round(y)]) {
                        depthBuffer[round(x)][round(y)] =  startDepth;
                        uint32_t Colour =
                                (255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
                       // if (isTriangleInsideViewFrustrum(triangle)) {
                        window.setPixelColour(round(x), round(y), Colour);
                       //       }else{break;}
                    }
                }
            }
            x1 += invSlope3;
            x2 += invSlope2;
            depth1 += depthSlope3;
            depth2 += depthSlope2;
        }


    drawRenderLine(triangle[0], triangle[1], LineColour, window, depthBuffer);
    drawRenderLine(triangle[1], triangle[2], LineColour, window, depthBuffer);
    drawRenderLine(triangle[2], triangle[0], LineColour, window, depthBuffer);
}






void RenderScene(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles,glm::vec3 cameraPosition,mat3 Camera_Orientation,std::vector<std::vector<float>> depthBuffer) {
    window.clearPixels();
    //int width = WIDTH;
   // int height = HEIGHT;
    //std::vector<std::vector<float>> depthBuffer(width, std::vector<float>(height, 0.001));
   // std::vector<std::vector<float>> depthBuffer(width, std::vector<float>(height, -(numeric_limits<float>::infinity())));


    //glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
    float focalLength = 2.0f;
    // Image plane scaling factor
    float scalingFactor = 150.0f;
    for (const ModelTriangle& modelTriangle : modelTriangles) {
        CanvasTriangle canvasTriangle;
        for (int i = 0; i < 3; i++) {
            vec3 vertexPosition = modelTriangle.vertices[i];
            glm::vec3 canvasPoint = getCanvasIntersectionPoint(cameraPosition, vertexPosition, focalLength, scalingFactor,Camera_Orientation);
           // if (canvasPoint.x != -1 && canvasPoint.y != -1) {
                CanvasPoint p = {canvasPoint.x, canvasPoint.y,canvasPoint.z};
                canvasTriangle.vertices[i] = p;
           // }
        }
        //drawSpecialTriangle(canvasTriangle,modelTriangle.colour,window);
        //drawFilledTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window);
        //if (isTriangleInsideViewFrustrum(canvasTriangle)) {
            drawRenderTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window, depthBuffer);
        //}
    }
}








