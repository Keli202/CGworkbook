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

using namespace std;
using namespace glm;



mat3 Camera_Orientation(vec3( -1.0,    0.0,    0.0),vec3( 0.0, 1.0,0.0),vec3( 0.0, 0.0, 1.0));
vec3 center = vec3(0.0f, 0.0f, 0.0f);
mat3 rotation_y(float t) {return mat3(vec3( cos(t), 0.0, sin(t)),vec3(    0.0, 1.0,    0.0),vec3(-sin(t), 0.0, cos(t)));}
mat3 rotation_x(float t) {return mat3(vec3( 1.0,    0.0,    0.0),vec3( 0.0, cos(t),-sin(t)),vec3( 0.0, sin(t), cos(t)));}
vec3 lightPosition(0.0f, 0.6f, 1.0f);
std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, std::numeric_limits<float>::infinity()));







void ClearDepthBuffer() {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            depthBuffer[x][y] = std::numeric_limits<float>::infinity();
        }
    }
}



vec3 getCanvasIntersectionPoint(vec3 cameraPosition,vec3 vertexPosition,float focalLength, float scalingFactor,mat3 Camera_Orientation){
    glm::vec3 relativePosition = vertexPosition - cameraPosition;
    relativePosition =  Camera_Orientation*relativePosition;
    if (-relativePosition.z <0) {
        return vec3(-1, -1, -1);// Vertex is behind the camera

    }
    // Calculate the projection coordinates of the vertices on the image plane
    double u = (focalLength * relativePosition.x / relativePosition.z) * scalingFactor + (WIDTH / 2);
    double v = (focalLength * relativePosition.y / relativePosition.z) * scalingFactor + (HEIGHT / 2);
    float depth = 1.0f / (relativePosition.z);
    return vec3(u, v, depth);
}

void drawRenderLine(CanvasPoint from, CanvasPoint to, Colour inputColour, DrawingWindow &window,vector<vector<float>>& depthBuffer){
    float xdistance = to.x - from.x;
    float ydistance = to.y - from.y;
    float zdistance = to.depth-from.depth;

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

        if (roundedX < 0 || roundedX >= WIDTH || roundedY < 0 || roundedY >= HEIGHT) {
            continue;
        }
        if (roundedX > 0 && roundedX <=WIDTH && roundedY > 0 && roundedY < HEIGHT) {
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

                if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT && depth < depthBuffer[x][y]) {
                    depthBuffer[x][y] = depth;
                    //float u = lambda0 * v0.texturePoint.y + lambda1 * v1.texturePoint.x + lambda2 * v2.texturePoint.u;
                    //float v = lambda0 * v0.texturePoint.y + lambda1 * v1.texturePoint.x + lambda2 * v2.texturePoint.v;

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
bool wire = false;

void RenderScene(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles,glm::vec3 cameraPosition,mat3 Camera_Orientation,std::vector<std::vector<float>> depthBuffer) {
    window.clearPixels();
    float focalLength = 2.0f;
    float scalingFactor = 150.0f;
    for (const ModelTriangle& modelTriangle : modelTriangles) {
        CanvasTriangle canvasTriangle;
        bool isValidTriangle = true;
        for (int i = 0; i < 3; i++) {
            vec3 vertexPosition = modelTriangle.vertices[i];
            glm::vec3 canvasPoint = getCanvasIntersectionPoint(cameraPosition, vertexPosition, focalLength, scalingFactor,Camera_Orientation);
            if (canvasPoint.x == -1 && canvasPoint.y == -1) {
                isValidTriangle = false;
                break;
            }
            CanvasPoint p = {canvasPoint.x, canvasPoint.y,canvasPoint.z};
            canvasTriangle.vertices[i] = p;
        }
        if (isValidTriangle) {
            if(wire){
                drawRenderLine(canvasTriangle.v0(), canvasTriangle.v1(), modelTriangle.colour, window, depthBuffer);
                drawRenderLine(canvasTriangle.v1(), canvasTriangle.v2(), modelTriangle.colour, window, depthBuffer);
                drawRenderLine(canvasTriangle.v2(), canvasTriangle.v0(), modelTriangle.colour, window, depthBuffer);
            }else{
                drawRenderTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window, depthBuffer);
            }
        }

    }
}


void readOBJColour(const std::string& filename, std::vector<std::pair<std::string, Colour>>& palette) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open material file " << filename << std::endl;
        return;
    }
    std::string line;
    std::string currentMaterial;

    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "newmtl") {
            currentMaterial = tokens[1];
        } else if (tokens[0] == "Kd") {
            int red = static_cast<int>(std::stof(tokens[1]) * 255);
            int green = static_cast<int>(std::stof(tokens[2]) * 255);
            int blue = static_cast<int>(std::stof(tokens[3]) * 255);
            palette.emplace_back(currentMaterial, Colour(red, green, blue));
        } else if (tokens[0] == "map_Kd") {
            palette.back().second.name = tokens[1];
        }
    }
    file.close();
}

std::vector<ModelTriangle> readOBJ(const std::string& filename, float scale) {
    std::ifstream file(filename);
    std::vector<ModelTriangle> triangles;
    std::vector<vec3> vertices;
    std::vector<TexturePoint> texturePoints;
    std::vector<std::pair<std::string, Colour>> palette;
    std::string textureName;

    readOBJColour("../textured-cornell-box.mtl", palette);

    std::string currentMaterial;
    std::string line;

    while (getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "v") {
            vec3 vertex(std::stof(tokens[1]) * scale, std::stof(tokens[2]) * scale, std::stof(tokens[3]) * scale);
            vertices.push_back(vertex);
        } else if (tokens[0] == "vt") {
            TexturePoint texPoint(std::stof(tokens[1]), std::stof(tokens[2]));
            texturePoints.push_back(texPoint);
        } else if (tokens[0] == "f") {
            std::vector<std::string> v1 = split(tokens[1], '/');
            std::vector<std::string> v2 = split(tokens[2], '/');
            std::vector<std::string> v3 = split(tokens[3], '/');
            //Colour triangleColour = palette.back().second;
            Colour triangleColour(255, 255, 255);
            for (const auto& item : palette) {
                if (item.first == currentMaterial) {
                    triangleColour = item.second;
                    textureName = item.second.name;
                    break;
                }
            }
            ModelTriangle triangle(
                    vertices[stoi(v1[0]) - 1],
                    vertices[stoi(v2[0]) - 1],
                    vertices[stoi(v3[0]) - 1],
                    triangleColour
            );
            if (!v1[1].empty() && !v2[1].empty() && !v3[1].empty() && v1.size() > 1 && v2.size() > 1 && v3.size() > 1) {

                triangle.texturePoints[0] = texturePoints[stoi(v1[1]) - 1];
                triangle.texturePoints[1] = texturePoints[stoi(v2[1]) - 1];
                triangle.texturePoints[2] = texturePoints[stoi(v3[1]) - 1];
            }
            //triangle.normal= calculateNormal(triangle);
            triangles.push_back(triangle);
        } else if (tokens[0] == "usemtl") {
            currentMaterial = tokens[1];
        }
    }

    file.close();
    return triangles;
}










void handleEvent(const std::vector<ModelTriangle>& modelTriangles,vec3& cameraPosition,SDL_Event event, DrawingWindow &window,std::vector<std::vector<float>> depthBuffer) {
    cout<<"1"<<endl;
    if (event.type == SDL_KEYDOWN) {
        cout<<"2"<<endl;

        if (event.key.keysym.sym == SDLK_1) {cout<<"3"<<endl;
            wire=!wire;}
        else if (event.key.keysym.sym == SDLK_2) { RenderScene(window,modelTriangles,cameraPosition,Camera_Orientation,depthBuffer);}


    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    window.clearPixels();
    if(wire){RenderScene(window,modelTriangles,cameraPosition,Camera_Orientation,depthBuffer);}

}













int main(int argc, char *argv[]) {


    glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
    std::vector<ModelTriangle> modelTriangles=readOBJ("../cornell-box.obj", 0.35);
    //std::vector<ModelTriangle> modelTriangles=readOBJ("../textured-cornell-box.obj", 0.35);
    //std::vector<ModelTriangle> modelTriangles=readOBJ("../sphere.obj", 0.35);



    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;
    while (true) {
        window.clearPixels();
        ClearDepthBuffer();
        handleEvent(modelTriangles,cameraPosition,event,window,depthBuffer);

//cout<<"1"<<endl;
        if(window.pollForInputEvents(event)) {
            //handleEvent(modelTriangles,cameraPosition,event,window,depthBuffer);
        }
//        if (orbitEnabled) {
//            float angle = -M_PI / 180;
//            orbitAndLookAt(cameraPosition, center, angle);
//        }
        window.renderFrame();
    }


}