
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

vec3 getCanvasIntersectionPoint(vec3 cameraPosition,vec3 vertexPosition,float focalLength, float scalingFactor,mat3 Camera_Orientation){

    glm::vec3 relativePosition = vertexPosition - cameraPosition;
    //relativePosition =  Camera_Orientation*relativePosition+cameraPosition;
    relativePosition =  Camera_Orientation*relativePosition;
    if (-relativePosition.z<=0) {
        return vec3(-1, -1, -1);
    }
    // Calculate the projection coordinates of the vertices on the image plane
    double u = (focalLength * relativePosition.x / relativePosition.z)*scalingFactor + (WIDTH / 2);
    //cout<<"X0:"<<relativePosition.x<<endl;
    //cout<<"y0:"<<relativePosition.y<<endl;
    double v = (focalLength * relativePosition.y / relativePosition.z)*scalingFactor + (HEIGHT / 2);
    if (u >= 0 && u < WIDTH && v >= 0 && v < HEIGHT) {
        float depth = 1.0f/(relativePosition.z);
        return vec3(u, v,depth);
    }
    else {
        if (u < 0) u = 0;
        if (u >= WIDTH) u = WIDTH - 1;
        if (v < 0) v = 0;
        if (v >= HEIGHT) v = HEIGHT - 1;
        float depth = 1.0f / (relativePosition.z);
        return vec3(u, v, depth);
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

        int roundedX = round(x);
        int roundedY = round(y);

//        if (roundedX < 0 || roundedX >= window.width || roundedY < 0 || roundedY >= window.height) {
//            continue;
//        }
        if (roundedX > 0 && roundedX < WIDTH && roundedY > 0 && roundedY < HEIGHT) {

            if (z < depthBuffer[round(x)][round(y)]) {
                uint32_t Colour = (255 << 24) + (inputColour.red << 16) + (inputColour.green << 8) + inputColour.blue;
                window.setPixelColour(round(x), round(y), Colour);
                depthBuffer[round(x)][round(y)] = z;
            }
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
    int minX = min({v0.x, v1.x, v2.x});
    int maxX = max({v0.x, v1.x, v2.x});
    int minY = min({v0.y, v1.y, v2.y});
    int maxY = max({v0.y, v1.y, v2.y});

    // Iterate over the bounding box
        for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) {
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
                        uint32_t colour =
                                (255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
                        window.setPixelColour(x, y, colour);

                }
            }
        }
    }

    drawRenderLine(v0, v1, LineColour, window, depthBuffer);
    drawRenderLine(v1, v2, LineColour, window, depthBuffer);
    drawRenderLine(v2, v0, LineColour, window, depthBuffer);
}


void RenderScene(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles,glm::vec3 cameraPosition,mat3 Camera_Orientation,std::vector<std::vector<float>> depthBuffer) {
    window.clearPixels();

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
            drawRenderTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window, depthBuffer);

    }
}








