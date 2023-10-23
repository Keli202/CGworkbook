//

//
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
#ifndef REDNOISE_DRAWTRIANGLE_H
#define REDNOISE_DRAWTRIANGLE_H

#endif //REDNOISE_DRAWTRIANGLE_H
//Week2
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

//Week2
void drawTriangle(SDL_Event event ,Colour colour,DrawingWindow &window){

    // Generate random triangle vertices and colors
    CanvasPoint vertex1(rand() % window.width, rand() % window.height);
    CanvasPoint vertex2(rand() % window.width, rand() % window.height);
    CanvasPoint vertex3(rand() % window.width, rand() % window.height);

    // Create the CanvasTriangle object and draw the triangle
    CanvasTriangle triangle({vertex1, vertex2, vertex3});
    drawLine(triangle[0], triangle[1], colour, window);
    drawLine(triangle[1], triangle[2], colour, window);
    drawLine(triangle[2], triangle[0], colour, window);

}
void drawKeyTriangle(SDL_Event event ,Colour colour,DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_u) {
            drawTriangle(event, colour, window);
        } else if (event.type == SDL_MOUSEBUTTONDOWN) {
            window.savePPM("output.ppm");
            window.saveBMP("output.bmp");
        }
    }
}

//week2
// Function to draw a filled triangle with specified color and stroked triangle over it
void drawFilledTriangles(CanvasTriangle triangle, Colour fillColour, Colour LineColour, DrawingWindow &window) {
    // Sorting vertices by y-coordinate
    if (triangle[0].y > triangle[1].y) std::swap(triangle[0], triangle[1]);
    if (triangle[0].y > triangle[2].y) std::swap(triangle[0], triangle[2]);
    if (triangle[1].y > triangle[2].y) std::swap(triangle[1], triangle[2]);

    // Calculate slopes of the three edges
    float slope1 = (triangle[1].x - triangle[0].x) / (triangle[1].y - triangle[0].y);
    float slope2 = (triangle[2].x - triangle[0].x) / (triangle[2].y - triangle[0].y);
    float slope3 = (triangle[2].x - triangle[1].x) / (triangle[2].y - triangle[1].y);

    // Initialize x coordinates for the three edges
    float x1 = triangle[0].x;
    float x2 = triangle[0].x;

    //uint32_t Colour = (255 << 24) + (fillColour.red << 16) + (fillColour.green << 8) + fillColour.blue;
    // Loop through each scanline (row) of the triangle
    for (float y = triangle[0].y; y < triangle[1].y; y++) {
        // Draw a horizontal line between x1 and x2 with fill colour
        CanvasPoint from(round(x1),round(y)) ;
        CanvasPoint to (round(x2),round(y)) ;
        drawLine(from,to,fillColour,window);
        // Update x coordinates for the three edges
        x1 += slope1;
        x2 += slope2;
    }
    float x3 = triangle[1].x;

    for (float y = triangle[1].y ; y <= triangle[2].y; y++) {
        CanvasPoint from(round(x3),round(y)) ;
        CanvasPoint to (round(x2),round(y)) ;
        drawLine(from,to,fillColour,window);
        x3 += slope3;
        x2 += slope2;
    }
    // Draw the stroked triangle (outline)
    drawLine(triangle[0], triangle[1], LineColour, window);
    drawLine(triangle[1], triangle[2], LineColour, window);
    drawLine(triangle[2], triangle[0], LineColour, window);
}
//week2
void drawSpecialTriangle(CanvasTriangle triangle,Colour colour,DrawingWindow &window){
    drawLine(triangle[0], triangle[1], colour, window);
    drawLine(triangle[1], triangle[2], colour, window);
    drawLine(triangle[2], triangle[0], colour, window);

}

//week2
void drawTexturedTriangle(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window) {
    CanvasPoint top = triangle.v0();
    CanvasPoint middle = triangle.v1();
    CanvasPoint bottom = triangle.v2();
    // Sort vertices by y-coordinate
    if (top.y > middle.y) swap(top, middle);
    if (top.y > bottom.y) swap(top, bottom);
    if (middle.y > bottom.y) swap(middle, bottom);
    // Calculate slopes for the top and bottom edges
    float slope1 = (middle.x - top.x) / (middle.y - top.y);
    float slope2 = (bottom.x - top.x) / (bottom.y - top.y);

    float left_x = top.x;
    float right_x = top.x;
    // First part of triangle (top to middle)
    for (int y = top.y; y <= middle.y; ++y) {
        // Calculate the corresponding texturePoints for CanvasPoints
        //Calculate the ratio between the canvas and texture
        //Left ratio
        float ratio1 = (y - top.y) / (middle.y - top.y);
        //Right ratio
        float ratio2 = (y - top.y) / (bottom.y - top.y);
        float texXLeft = top.texturePoint.x + (middle.texturePoint.x - top.texturePoint.x) * ratio1;
        float texXRight = top.texturePoint.x + (bottom.texturePoint.x - top.texturePoint.x) * ratio2;
        float texYLeft = top.texturePoint.y + (middle.texturePoint.y - top.texturePoint.y) * ratio1;
        float texYRight = top.texturePoint.y + (bottom.texturePoint.y - top.texturePoint.y) * ratio2;

        // Draw the horizontal line with texture mapping
        for (int x = round(left_x); x <= round(right_x); ++x) {
            // Interpolate t based on the current x position
            //Calculate the ratio between the horizontal coordinates
            float ratio = (x - left_x) / (right_x - left_x);
            int texX = texXLeft + (texXRight - texXLeft) * ratio;
            int texY = texYLeft + (texYRight - texYLeft) * ratio;

            //Ensure that texture coordinates are within valid ranges to avoid out-of-bounds errors when accessing texture data.
            texX = std::max(0, std::min(static_cast<int>(texX), static_cast<int>(texture.width) - 1));
            texY = std::max(0, std::min(static_cast<int>(texY), static_cast<int>(texture.height) - 1));
            // Calculate the index to access the texture pixel
            int index = texY * texture.width + texX;
            window.setPixelColour(x, y, texture.pixels[index]);
        }
        left_x += slope1;
        right_x += slope2;
    }

    // Initialize left and right x-coordinates for the bottom edge
    left_x = middle.x;
    float slope3 = (bottom.x - middle.x) / (bottom.y - middle.y);

    // Second part of triangle (middle to bottom)
    for (int y = middle.y + 1; y <= bottom.y; ++y) {
        // Calculate the corresponding texture coordinates for Canvas Points
        float ratio3 = (y - middle.y) / (bottom.y - middle.y);
        float ratio2 = (y - top.y) / (bottom.y - top.y);
        float texXLeft2 = middle.texturePoint.x + (bottom.texturePoint.x - middle.texturePoint.x) * ratio3;
        float texXRight = top.texturePoint.x + (bottom.texturePoint.x - top.texturePoint.x) * ratio2;
        float texYLeft2 = middle.texturePoint.y + (bottom.texturePoint.y - middle.texturePoint.y) * ratio3;
        float texYRight = top.texturePoint.y + (bottom.texturePoint.y - top.texturePoint.y) * ratio2;

        // Draw the horizontal line with texture mapping
        for (int x = round(left_x); x <= round(right_x); ++x) {
            // Interpolate t based on the current x position
            float ratio = (x - left_x) / (right_x - left_x);
            int texX = texXLeft2 + (texXRight - texXLeft2) * ratio;
            int texY = texYLeft2 + (texYRight - texYLeft2) * ratio;
            texX = std::max(0, std::min(static_cast<int>(texX), static_cast<int>(texture.width) - 1));
            texY = std::max(0, std::min(static_cast<int>(texY), static_cast<int>(texture.height) - 1));
            // Calculate the index to access the texture pixel
            int index = texY * texture.width + texX;
            window.setPixelColour(x, y, texture.pixels[index]);
        }
        left_x += slope3;
        right_x += slope2;
    }
    drawSpecialTriangle(triangle,Colour(255,255,255),window);
}
