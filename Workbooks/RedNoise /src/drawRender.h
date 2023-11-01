
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
float interpolate(float from, float to, float factor) {
    return from + (to - from) * factor;
}

void drawRenderLine(CanvasPoint from, CanvasPoint to, Colour inputColour, DrawingWindow &window,vector<vector<float>>& depthBuffer){
    float xdistance = to.x - from.x;
    float ydistance = to.y - from.y;
    float zdistance = to.depth-from.depth;
    //float depthFrom = 1.0 / from.depth;
    //float depthTo = 1.0 / to.depth;
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
            if (z<depthBuffer[round(x)][round(y)]){
           // cout<<"dbHHHHH:"<<depthBuffer[round(x)][round(y)]<<endl;
            //cout<<"zHHHH:"<<z<<endl;
            uint32_t Colour = (255 << 24) + (inputColour.red << 16) + (inputColour.green << 8) + inputColour.blue;
            window.setPixelColour(round(x), round(y), Colour);
            depthBuffer[round(x)][round(y)]=z;
            // cout<<"db2:"<<depthBuffer[round(x)][round(y)]<<endl;
            // cout<<"z2:"<<z<<endl;
        }
    }
}

void drawRenderTriangles(CanvasTriangle triangle, Colour fillColour, Colour LineColour, DrawingWindow &window,vector<vector<float>>& depthBuffer) {
    // Sort the vertices by y-coordinate
    if (triangle[0].y > triangle[1].y) std::swap(triangle[0], triangle[1]);
    if (triangle[0].y > triangle[2].y) std::swap(triangle[0], triangle[2]);
    if (triangle[1].y > triangle[2].y) std::swap(triangle[1], triangle[2]);

    CanvasPoint v0 = triangle[0];
    CanvasPoint v1 = triangle[1];
    CanvasPoint v2 = triangle[2];

    // Compute bounding box of the triangle
    int minX = std::min({v0.x, v1.x, v2.x});
    int maxX = std::max({v0.x, v1.x, v2.x});
    int minY = std::min({v0.y, v1.y, v2.y});
    int maxY = std::max({v0.y, v1.y, v2.y});

    // Iterate over the bounding box
    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            if (x < 0 || x >= window.width || y < 0 || y >= window.height) {
                continue;
            }
            // Compute Barycentric coordinates
            float lambda0 = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) /
                            ((v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y));
            float lambda1 = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) /
                            ((v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y));
            float lambda2 = 1.0f - lambda0 - lambda1;

            // Check if the point is inside the triangle
            if (lambda0 >= 0 && lambda1 >= 0 && lambda2 >= 0) {
                // Interpolate depth
                float depth = lambda0 * v0.depth + lambda1 * v1.depth + lambda2 * v2.depth;
                // Depth test
                if (depth < depthBuffer[x][y]) {
                    depthBuffer[x][y] = depth;
                    uint32_t colour = (255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
                    window.setPixelColour(x, y, colour);
                }
            }
        }
    }

    drawRenderLine(v0, v1, LineColour, window, depthBuffer);
    drawRenderLine(v1, v2, LineColour, window, depthBuffer);
    drawRenderLine(v2, v0, LineColour, window, depthBuffer);
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
//        //leftDepth = leftDepth2 + (y - middle.y) * depthSlope3;
//        //rightDepth = rightDepth2 + (y - middle.y) * depthSlope2;
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


//void drawRenderTriangles(CanvasTriangle triangle, Colour fillColour, Colour LineColour, DrawingWindow &window, vector<vector<float>>& depthBuffer) {
//
//    if (triangle[0].y > triangle[1].y) std::swap(triangle[0], triangle[1]);
//    if (triangle[0].y > triangle[2].y) std::swap(triangle[0], triangle[2]);
//    if (triangle[1].y > triangle[2].y) std::swap(triangle[1], triangle[2]);
//
//    float invSlope1 = (triangle[1].x - triangle[0].x) / (triangle[1].y - triangle[0].y);
//    float invSlope2 = (triangle[2].x - triangle[0].x) / (triangle[2].y - triangle[0].y);
//    float depthSlope1 = (triangle[1].depth - triangle[0].depth) / (triangle[1].y - triangle[0].y);
//    float depthSlope2 = (triangle[2].depth - triangle[0].depth) / (triangle[2].y - triangle[0].y);
//
//    float x1 = triangle[0].x;
//    float x2 = triangle[0].x;
//    float depth1 = triangle[0].depth;
//    float depth2 = triangle[0].depth;
//        for (float y = triangle[0].y; y < triangle[1].y; y++) {
//            float startX = std::min(x1, x2);
//            float endX = std::max(x1, x2);
//            float startDepth = std::min(depth1, depth2);
//            float depthFrom = depth1;
//            float depthTo =  depth2;
//            for (float x = round(startX); x <= round(endX); x++) {
//                if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
//                    float t = abs((x - startX) / (endX - startX));
//                    if (endX == startX) {
//                        t = 0;
//                    } else {
//                        t = clamp(t, 0.0f, 1.0f);
//                    }
//                    //cout<<"t1:"<<t<<endl;
//                    float currentDepth = interpolate(depthFrom, depthTo, t);
//                    //if (-1.0f / startDepth > depthBuffer[round(x)][round(y)]) {
//                        if (currentDepth <depthBuffer[round(x)][round(y)]) {
//                        depthBuffer[round(x)][round(y)] = currentDepth;
//                        uint32_t Colour =
//                                (255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
//                       // if (isTriangleInsideViewFrustrum(triangle)) {
//                        window.setPixelColour(round(x), round(y), Colour);
//                        //  }else{break;}
//                    }
//                }
//            }
//            x1 += invSlope1;
//            x2 += invSlope2;
//            depth1 += depthSlope1;
//            depth2 += depthSlope2;
//        }
//
//        float invSlope3 = (triangle[2].x - triangle[1].x) / (triangle[2].y - triangle[1].y);
//        float depthSlope3 = (triangle[2].depth - triangle[1].depth) / (triangle[2].y - triangle[1].y);
//
//        x1 = triangle[1].x;
//        depth1 = triangle[1].depth;
//
//        for (float y = triangle[1].y; y <= triangle[2].y; y++) {
//            float startX = std::min(x1, x2);
//            float endX = std::max(x1, x2);
//            float startDepth = std::min(depth1, depth2);
//            float depthFrom =  depth1;
//            float depthTo =  depth2;
//            for (float x = round(startX); x <= round(endX); x++) {
//                if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
//                    float t = abs((x - startX) / (endX - startX));
//                    if (endX == startX) {
//                        t = 0;
//                    } else {
//                        t = clamp(t, 0.0f, 1.0f);
//                    }
//                   // cout<<"t2:"<<t<<endl;
//                    float currentDepth = interpolate(depthFrom, depthTo, t);
//                    //if (-1.0f / startDepth > depthBuffer[round(x)][round(y)]) {
//                        if ( currentDepth < depthBuffer[round(x)][round(y)]) {
//                        depthBuffer[round(x)][round(y)] =  currentDepth;
//                        uint32_t Colour =(255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
//                       // if (isTriangleInsideViewFrustrum(triangle)) {
//                        window.setPixelColour(round(x), round(y), Colour);
//                       //       }else{break;}
//                    }
//                }
//            }
//            x1 += invSlope3;
//            x2 += invSlope2;
//            depth1 += depthSlope3;
//            depth2 += depthSlope2;
//        }
//
//
//    drawRenderLine(triangle[0], triangle[1], LineColour, window, depthBuffer);
//    drawRenderLine(triangle[1], triangle[2], LineColour, window, depthBuffer);
//    drawRenderLine(triangle[2], triangle[0], LineColour, window, depthBuffer);
//}






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








