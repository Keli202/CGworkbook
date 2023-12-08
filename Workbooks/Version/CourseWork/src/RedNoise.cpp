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
#include <RayTriangleIntersection.h>
#include "RayTrace.h"
using namespace std;
using namespace glm;



mat3 Camera_Orientation(vec3( 1.0,    0.0,    0.0),vec3( 0.0, 1.0,0.0),vec3( 0.0, 0.0, 1.0));
vec3 center = vec3(0.0f, 0.0f, 0.0f);
mat3 rotation_y(float t) {return mat3(vec3( cos(t), 0.0, sin(t)),vec3(    0.0, 1.0,    0.0),vec3(-sin(t), 0.0, cos(t)));}
mat3 rotation_x(float t) {return mat3(vec3( 1.0,    0.0,    0.0),vec3( 0.0, cos(t),-sin(t)),vec3( 0.0, sin(t), cos(t)));}
//vec3 lightPosition(0.0f, 0.6f, 1.0f);
std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, std::numeric_limits<float>::infinity()));
enum RenderModel{
    None,
    WireFrame,
    Render,
    Diffuse,
    Texture,
    Specular,
    Soft,
    Phong,
    Mirror
};
TextureMap normalTexture("../texture.ppm");
static bool orbitEnabled = false;

//vector <vec3> lightSources() {
//    vec3 lightPosition(0.0f, 0.6f, 1.0f);
//    int numSamples = 10;
//    float sampleRadius = 0.1f;
//
//    for (int i = 0; i < numSamples; i++) {
//        float angle = (2 * M_PI / numSamples) * i;
//        lightSources().push_back(lightPosition + glm::vec3(sampleRadius * cos(angle), 0, sampleRadius * sin(angle)));
//    }
//    return lightSources();
//}

void moveLights(float scale, string axis){
    vector<vec3>lightSource=lightSources();
    for(int i=0;i < lightSource.size();i++){
        if(axis=="x"){
            lightSource[i].x+=scale;
            //cout<<"lights:"<<lightSource[i].x<<endl;
        }else if(axis=="y"){
            lightSource[i].y+=scale;
        }else if(axis=="z"){
            lightSource[i].z+=scale;
        }

    }
}

void ClearDepthBuffer() {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            depthBuffer[x][y] = std::numeric_limits<float>::infinity();
        }
    }
}



vec3 getCanvasIntersectionPoint(vec3 cameraPosition,vec3 vertexPosition,float focalLength, float scalingFactor,mat3 Camera_Orientation){
    glm::vec3 relativePosition = vertexPosition - cameraPosition;
    //relativePosition =  Camera_Orientation*relativePosition;
    relativePosition =  relativePosition*Camera_Orientation;

    if (-relativePosition.z <0) {
        return vec3(-1, -1, -1);// Vertex is behind the camera

    }
    // Calculate the projection coordinates of the vertices on the image plane
    double u = -(focalLength * relativePosition.x / relativePosition.z) * scalingFactor + (WIDTH / 2);
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
        if (roundedX >= 0 && roundedX <WIDTH && roundedY >= 0 && roundedY < HEIGHT) {
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




void drawTextureTriangle(CanvasTriangle triangle, TextureMap& texture, DrawingWindow &window, vector<vector<float>>& depthBuffer) {
    if (triangle[0].y > triangle[1].y) std::swap(triangle[0], triangle[1]);
    if (triangle[0].y > triangle[2].y) std::swap(triangle[0], triangle[2]);
    if (triangle[1].y > triangle[2].y) std::swap(triangle[1], triangle[2]);
    CanvasPoint v0 = triangle[0];
    CanvasPoint v1 = triangle[1];
    CanvasPoint v2 = triangle[2];
    int minX = min({v0.x, v1.x, v2.x});
    int maxX = max({v0.x, v1.x, v2.x});
    int minY = min({v0.y, v1.y, v2.y});
    int maxY = max({v0.y, v1.y, v2.y});
    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            if (x < 0 || x >= window.width || y < 0 || y >= window.height) continue;
            float lambda0 = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) /
                            ((v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y));
            float lambda1 = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) /
                            ((v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y));
            float lambda2 = 1.0f - lambda0 - lambda1;
            if (lambda0 >= 0 && lambda1 >= 0 && lambda2 >= 0) {
                float depth = lambda0 * v0.depth + lambda1 * v1.depth + lambda2 * v2.depth;
                if (depth < depthBuffer[x][y]) {
                    depthBuffer[x][y] = depth;
                    float texX = lambda0 * v0.texturePoint.x + lambda1 * v1.texturePoint.x + lambda2 * v2.texturePoint.x;
                    float texY = lambda0 * v0.texturePoint.y + lambda1 * v1.texturePoint.y + lambda2 * v2.texturePoint.y;
                    int texXInt = static_cast<int>(std::round(texX)) % texture.width;
                    int texYInt = static_cast<int>(std::round(texY)) % texture.height;
                    int textureIndex = texYInt * texture.width + texXInt;
                    uint32_t colour = texture.pixels[textureIndex];
                    window.setPixelColour(x, y, colour);
                }
            }
        }
    }
}





void RenderScene(DrawingWindow &window, const std::vector<ModelTriangle>& modelTriangles,glm::vec3 cameraPosition,mat3 Camera_Orientation,std::vector<std::vector<float>> depthBuffer,bool ifwire) {
    window.clearPixels();
    float focalLength = 2.0f;
    float scalingFactor = 150.0f;
    for (const ModelTriangle& modelTriangle : modelTriangles) {
        CanvasTriangle canvasTriangle;
        bool isValidTriangle = true;
        for (int i = 0; i < 3; i++) {
            vec3 vertexPosition = modelTriangle.vertices[i];
            glm::vec3 canvasPoint = getCanvasIntersectionPoint(cameraPosition, vertexPosition, focalLength, scalingFactor,Camera_Orientation);
            if (canvasPoint.x <0 && canvasPoint.y <0) {
                isValidTriangle = false;
                break;
            }
            CanvasPoint p = {canvasPoint.x, canvasPoint.y,canvasPoint.z};
            canvasTriangle.vertices[i] = p;
            canvasTriangle.vertices[i].texturePoint=modelTriangle.texturePoints[i];

        }
        if (isValidTriangle) {
            if(ifwire){
                drawRenderLine(canvasTriangle.v0(), canvasTriangle.v1(), modelTriangle.colour, window, depthBuffer);
                drawRenderLine(canvasTriangle.v1(), canvasTriangle.v2(), modelTriangle.colour, window, depthBuffer);
                drawRenderLine(canvasTriangle.v2(), canvasTriangle.v0(), modelTriangle.colour, window, depthBuffer);
            }else{
                if (!modelTriangle.colour.name.empty()) {
                    TextureMap texture=normalTexture;
                    for(int j = 0; j < canvasTriangle.vertices.size(); j++) {
                        canvasTriangle.vertices[j].texturePoint.x *= texture.width;
                        canvasTriangle.vertices[j].texturePoint.y *= texture.height;
                    }
                    drawTextureTriangle(canvasTriangle, texture, window, depthBuffer);
                }else {
                    drawRenderTriangles(canvasTriangle, modelTriangle.colour, modelTriangle.colour, window,
                                        depthBuffer);
                }
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
            float r = std::stof(tokens[1]);
            float g = std::stof(tokens[2]);
            float b = std::stof(tokens[3]);
            int red = static_cast<int>(r * 255);
            int green = static_cast<int>(g * 255);
            int blue = static_cast<int>(b * 255);
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

//    readOBJColour("../textured-cornell-box.mtl", palette);
    readOBJColour("../cornell-box.mtl", palette);

    std::string currentMaterial;
    std::string line;
    bool mirror=false;
    bool glass=false;
    bool metal=false;
    bool texture=false;

    string material;
    while (getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "v") {
            vec3 vertex(std::stof(tokens[1]) * scale, std::stof(tokens[2]) * scale, std::stof(tokens[3]) * scale);
            vertices.push_back(vertex);
        } else if (tokens[0]=="o"){
            material=tokens[1];
        }else if (tokens[0] == "vt") {
            TexturePoint texPoint(std::stof(tokens[1]), std::stof(tokens[2]));
            texturePoints.push_back(texPoint);
        } else if (tokens[0] == "f") {
            std::vector<std::string> v1 = split(tokens[1], '/');
            std::vector<std::string> v2 = split(tokens[2], '/');
            std::vector<std::string> v3 = split(tokens[3], '/');
            Colour triangleColour(255, 255, 255);
            for (const auto& item : palette) {
                if (item.first == currentMaterial) {
                    triangleColour = item.second;
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
            triangle.mirror = mirror;
            triangle.glass = glass;
            triangle.texture =texture;
//            cout<<"check"<<texture<<endl;
            triangle.metal=metal;
            //triangle.normal= calculateNormal(triangle);
            triangle.material=material;
            //cout<<"m:"<<triangle.material<<endl;
            triangles.push_back(triangle);
        } else if (tokens[0] == "usemtl") {
            mirror= (tokens[1]=="Mirror");
            glass=(tokens[1]=="Glass");
            metal=(tokens[1]=="Metal");
            texture=(tokens[1]=="Texture");
            //cout<<"t:"<<texture<<endl;
            currentMaterial = tokens[1];
            //cout<<"check"<<tokens[1]<<endl;

//            if(currentMaterial=="Texture"){
//                texture=true;
//                //cout<<"1"<<endl;
//            }else {
//                //cout<<"2"<<endl;
//                 texture =false;}
        }
    }

    file.close();
    return triangles;
}


mat3 lookAt(vec3 cameraPosition, vec3 lookPoint) {
    vec3 forward = normalize(cameraPosition-lookPoint);
    vec3 vertical = vec3(0.0f, 1.0f, 0.0f);
    vec3 right = normalize(cross(vertical, forward));
    vec3 up = normalize(cross(forward, right));
    return mat3(right, up, forward);
}

void orbitAndLookAt(vec3& cameraPosition, vec3 orbitCenter, float angle) {
    vec3 relativePosition = cameraPosition-orbitCenter;

    float cosA = cos(angle);
    float sinA = sin(angle);

    vec3 newPosition;
    newPosition.x = relativePosition.x * cosA + relativePosition.z * sinA;
    newPosition.y = relativePosition.y;
    newPosition.z = -relativePosition.x * sinA + relativePosition.z * cosA;
    cameraPosition = orbitCenter + newPosition;

    Camera_Orientation = lookAt(cameraPosition, orbitCenter);
}


















void handleEvent(const std::vector<ModelTriangle>& modelTriangles,vec3& cameraPosition,SDL_Event event, DrawingWindow &window,std::vector<std::vector<float>> depthBuffer, RenderModel& model) {
    float a=0.1f;
    float t=M_PI/180;
    float lightMoveStep = 0.1f;
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_1) {model = WireFrame;}
        else if (event.key.keysym.sym == SDLK_2) { model = Render;}
        else if (event.key.keysym.sym == SDLK_3) { model = Diffuse;}
        else if (event.key.keysym.sym == SDLK_4) { model = Texture;}
        else if (event.key.keysym.sym == SDLK_5) { model = Specular;}
        else if (event.key.keysym.sym == SDLK_6) { model = Soft;}
        else if (event.key.keysym.sym == SDLK_7) { model = Mirror;}


        else if (event.key.keysym.sym == SDLK_a) { cameraPosition.x += a; }
        else if (event.key.keysym.sym == SDLK_d) { cameraPosition.x -= a; }
        else if (event.key.keysym.sym == SDLK_w) { cameraPosition.y -= a; }
        else if (event.key.keysym.sym == SDLK_s) { cameraPosition.y += a; }
        else if (event.key.keysym.sym == SDLK_o) { cameraPosition.z += a; }
        else if (event.key.keysym.sym == SDLK_p) { cameraPosition.z -= a; }
        else if (event.key.keysym.sym == SDLK_k) { cameraPosition = cameraPosition * rotation_x(-t); }
        else if (event.key.keysym.sym == SDLK_i) { cameraPosition = cameraPosition * rotation_x(t); }
        else if (event.key.keysym.sym == SDLK_l) { cameraPosition = cameraPosition * rotation_y(-t); }
        else if (event.key.keysym.sym == SDLK_j) { cameraPosition = cameraPosition * rotation_y(t); }
        else if (event.key.keysym.sym == SDLK_LEFT) { Camera_Orientation = Camera_Orientation * rotation_y(t); }
        else if (event.key.keysym.sym == SDLK_RIGHT) { Camera_Orientation = Camera_Orientation * rotation_y(-t); }
        else if (event.key.keysym.sym == SDLK_UP) { Camera_Orientation = Camera_Orientation * rotation_x(t); }
        else if (event.key.keysym.sym == SDLK_DOWN) { Camera_Orientation = Camera_Orientation * rotation_x(-t); }
        else if (event.key.keysym.sym == SDLK_z) { lightPosition.x -= lightMoveStep; }
        else if (event.key.keysym.sym == SDLK_x) { lightPosition.x += lightMoveStep; }
        else if (event.key.keysym.sym == SDLK_c) { lightPosition.y += lightMoveStep; }
        else if (event.key.keysym.sym == SDLK_v) { lightPosition.y -= lightMoveStep; }
        else if (event.key.keysym.sym == SDLK_b) { lightPosition.z -= lightMoveStep; }
        else if (event.key.keysym.sym == SDLK_n) { lightPosition.z += lightMoveStep; }
        else if (event.key.keysym.sym == SDLK_r) {
            float r = sqrt(cameraPosition.x * cameraPosition.x + cameraPosition.z * cameraPosition.z);
            float theta = atan2(cameraPosition.z, cameraPosition.x);
            float deltaTheta = 0.05f;
            theta += deltaTheta;
            cameraPosition.x = r * cos(theta);
            cameraPosition.z = r * sin(theta);
        } else if (event.key.keysym.sym == SDLK_g) {
            orbitEnabled = !orbitEnabled;
            if (orbitEnabled) {
                orbitAndLookAt(cameraPosition, center, M_PI / 180);
            }
        }

    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    //window.clearPixels();

}













int main(int argc, char *argv[]) {


    glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
    std::vector<ModelTriangle> modelTriangles1=readOBJ("../cornell-box.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles2=readOBJ("../textured-cornell-box.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles3=readOBJ("../sphere.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles4=readOBJ("../loveBox3.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles5=readOBJ("../loveBox.obj", 0.35);



    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    RenderModel model = None;


    bool animationCompleted = false;

    SDL_Event event;
   // while (true) {
   while (!animationCompleted) {

            if(window.pollForInputEvents(event)) {
            handleEvent(modelTriangles1,cameraPosition,event,window,depthBuffer,model);
            handleEvent(modelTriangles2,cameraPosition,event,window,depthBuffer,model);
            handleEvent(modelTriangles3,cameraPosition,event,window,depthBuffer,model);

        }
        window.clearPixels();
        ClearDepthBuffer();

        if (orbitEnabled) {
            float angle = -M_PI / 180;
            orbitAndLookAt(cameraPosition, center, angle);
        }
        modelTriangles1= getVertexNormals(modelTriangles1);
        modelTriangles2= getVertexNormals(modelTriangles2);
        modelTriangles3= getVertexNormals(modelTriangles3);
        modelTriangles4= getVertexNormals(modelTriangles4);
        modelTriangles5= getVertexNormals(modelTriangles5);


        //vec3 lightPosition(0,0,0);
            if (model == WireFrame) {
                RenderScene(window, modelTriangles1, cameraPosition, Camera_Orientation, depthBuffer, true);
            }
            else if (model == Render) {
                RenderScene(window, modelTriangles2, cameraPosition, Camera_Orientation, depthBuffer, false);
            }
            else if (model == Diffuse) {
                DrawRay(modelTriangles4, window, cameraPosition, Camera_Orientation,  256.0f);
            }
            else if (model == Texture) {
                DrawRay(modelTriangles5, window, cameraPosition, Camera_Orientation,  256.0f);
            }
            else if (model == Specular) {
                DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation, 256.0f);
            }
            else if (model == Soft) {
                DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
            }
            else if (model == Mirror) {
                DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);
            }


//        for (int frameCount = 0; frameCount < 1000; frameCount++) {
//            window.clearPixels();
//            int n_zero = 5;
//            RenderScene(window, modelTriangles1, cameraPosition, Camera_Orientation, depthBuffer, true);
//            string name = string(n_zero - to_string(frameCount).length(), '0') + to_string(frameCount);
//            window.savePPM("../output1/"+name+".ppm");
//            cout << "saved " << frameCount << endl;
//
//                orbitEnabled = !orbitEnabled;
//                if (orbitEnabled) {
//                    orbitAndLookAt(cameraPosition, center, M_PI / 180);
//                }
//
//
//        }


//       for (int frameCount = 0; frameCount < 2; frameCount++) {
//           window.clearPixels();
//           int n_zero = 5;
//           DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);
//
//           //RenderScene(window, modelTriangles1, cameraPosition, Camera_Orientation, depthBuffer, false);
//           string name = string(n_zero - to_string(frameCount).length(), '0') + to_string(frameCount);
//           window.savePPM("../output3ray/"+name+".ppm");
//           cout << "saved " << frameCount << endl;
//
////           if(frameCount < 10) {
////               //DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
////               moveLights(-0.05f, "x");
////
////           }
////        else if(frameCount < 30) {
////               //DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
////               moveLights(0.05f, "x");
////
////           }
////        else if(frameCount < 40) {
////              // DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
////               moveLights(-0.05f, "y");
////
////           }else if(frameCount < 60) {
////              // DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
////               moveLights(-0.05f, "z");
////
////           }
//
//
//       }

//
//
//
//        for (int frameCount = 0; frameCount < 200; frameCount++) {
//
//           for (int frame = 0; frame < frameCount; ++frame) {
//               window.clearPixels();
//        //renderRasterising(mirrorTriangle2,focalLength,window,400.0f);
//        string frameNumber = string(n_zero - to_string(frame).length(), '0') + to_string(frame);
//        window.savePPM("../output14/" + frameNumber + ".ppm");
//        cout << "saved " << frame << endl;
//        float speed = -3.8 * M_PI / 180;
//        cameraPosition = cameraPosition * mat3(
//                vec3(cos(speed), 0.0, sin(speed)),
//                vec3(0.0, 1.0, 0.0),
//                vec3(-sin(speed), 0.0, cos(speed))
//        );
//        cameraOrientation = cameraOrientation * mat3(
//                vec3(cos(speed), 0.0, sin(speed)),
//                vec3(0.0, 1.0, 0.0),
//                vec3(-sin(speed), 0.0, cos(speed))
//        );
//        look_At();
//    }
//    for (int frame = 0; frame < frameCount; ++frame) {
//        window.clearPixels();
//        renderRasterising(triangles, focalLength, window, 400.0f);
//        string frameNumber = string(n_zero - to_string(frame).length(), '0') + to_string(frame);
//        window.savePPM("../output3/" + frameNumber + ".ppm");
//        cout << "saved " << frame << endl;
//        if(frame < 12) {
//            cameraPosition.x -= 0.06;
//        }
//        else if(frame < 36) {
//            cameraPosition.x += 0.06;
//        }
//        else if(frame < 48) {
//            cameraPosition.x -= 0.06;
//        }
//    }
//
//    for (int frame = 0; frame < frameCount; ++frame) {
//        window.clearPixels();
//        renderRasterising(triangles, focalLength, window, 400.0f);
//        string frameNumber = string(n_zero - to_string(frame).length(), '0') + to_string(frame);
//        window.savePPM("../output4/" + frameNumber + ".ppm");
//        cout << "saved " << frame << endl;
//        if(frame < 12) {
//            cameraPosition.y -= 0.06;
//        }
//        else if(frame < 36) {
//            cameraPosition.y += 0.06;
//        }
//        else if(frame < 48) {
//            cameraPosition.y -= 0.06;
//        }
//    }
//    for (int frame = 0; frame < frameCount; ++frame) {
//        window.clearPixels();
//        renderRasterising(triangles, focalLength, window, 400.0f);
//        string frameNumber = string(n_zero - to_string(frame).length(), '0') + to_string(frame);
//        window.savePPM("../output2/" + frameNumber + ".ppm");
//        cout << "saved " << frame << endl;
//        if(frame < 12) {
//            cameraPosition.z += 0.06;
//        }
//        else if(frame < 36) {
//            cameraPosition.z -= 0.06;
//        }
//        else if(frame < 48) {
//            cameraPosition.z -= 0.06;
//        }
//    }
//
//    for (int frame = 0; frame < frameCount; ++frame) {
//        window.clearPixels();
//        renderRayTracedScene(triangles2,window,focalLength,400.0f);
//        string frameNumber = string(n_zero - to_string(frame).length(), '0') + to_string(frame);
//        window.savePPM("../output8/" + frameNumber + ".ppm");
//        cout << "saved " << frame << endl;
//        cameraPosition.x -= 0.005f;
//        cameraPosition.z -= 0.02f;
//        cameraPosition.y += 0.004f;
//    }
















           // animationCompleted = true;

        window.renderFrame();



    }


}