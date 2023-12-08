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
#include <algorithm>
using namespace std;
using namespace glm;


std::vector<glm::vec3> lightSources() {
    glm::vec3 center(0.0f,0.6f,0.8f);
    float radius=0.1f;
    int numSources=1;
    float height=0.6f;
    std::vector<glm::vec3> lightSources;
    for (int i = 0; i < numSources; ++i) {
        //AI help me calculate the formulas
        float angle = 2.0f * M_PI * i / numSources;
        float x = center.x + radius * cos(angle);
        float z = center.z + radius * sin(angle);
        lightSources.push_back(glm::vec3(x, height, z));
    }
    return lightSources;


}
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
TextureMap normalTexture("texture.ppm");
static bool orbitEnabled = false;




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
    readOBJColour("cornell-box.mtl", palette);

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














vec3 lightPosition(0.0f,0.6f,1.2f);
//vector <vec3> lightSource{
//
//        vec3(0.2f, 0.6f, 1.0f),
//        vec3(-0.2f, 0.6f, 1.0f),
//        vec3(0.2, 0.6, 1.2),
//         vec3(-0.2, 0.6, 1.2),
//         vec3(0.0, 0.6, 1.0),
//         vec3(0.0, 0.6, 1.2),
//         vec3(0.0, 0.6, 0.8),
//         vec3(-0.2, 0.6, 0.8),
//         vec3(0.1, 0.6, -0.1),
//         vec3(0.1, 0.6, 0.0),
////         vec3(0.1, 0.8, 0.1),
////         vec3(-0.1, 0.9, -0.1),
////         vec3(-0.1, 0.9, 0.0),
////         vec3(-0.1, 0.9, 0.1),
////         vec3(0.0, 0.9, -0.1),
////         vec3(0.0, 0.9, 0.0),
////         vec3(0.0, 0.9, 0.1),
////         vec3(0.1, 0.9, -0.1),
////         vec3(0.1, 0.9, 0.0),
////         vec3(0.1, 0.9, 0.1),
////        vec3(-0.2, 1.0, 1.8),
////        vec3(-0.2, 1.0, 2.0),
////        vec3(-0.2, 1.0, 2.2),
////        vec3(0.0, 1.0, 1.8),
////        vec3(0.0, 1.0, 2.0),
////        vec3(0.0, 1.0, 2.2),
////        vec3(0.2, 1.0, 1.8),
////        vec3(0.2, 1.0, 2.0),
////        vec3(0.2, 1.0, 2.2)
//
//};

TextureMap TextureFile1("texture_2.ppm");



glm::vec3 calculateNormal(const ModelTriangle& triangle,RayTriangleIntersection &intersection) {
    glm::vec3 edge1 = triangle.vertices[1] - triangle.vertices[0];
    glm::vec3 edge2 = triangle.vertices[2] - triangle.vertices[0];
    glm::vec3 normal = glm::cross(edge1, edge2);
    intersection.intersectedTriangle.normal= glm::normalize(normal);
    return intersection.intersectedTriangle.normal;
}

glm::vec3 calculateFaceNormal(const ModelTriangle& triangle) {
    glm::vec3 edge1 = triangle.vertices[1] - triangle.vertices[0];
    glm::vec3 edge2 = triangle.vertices[2] - triangle.vertices[0];
    glm::vec3 normal = glm::cross(edge1, edge2);
    return glm::normalize(normal);
}

//AI
vec3 refract(const vec3& incident, const vec3& normal, float indexOfRefraction) {
    float cosi = clamp(-1.0f, 1.0f, dot(incident, normal));
    float etai = 1, etat = indexOfRefraction;
    vec3 n = normal;
    if (cosi < 0) {
        cosi = -cosi;
    } else {
        std::swap(etai, etat);
        n= -normal;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(0,0,0) : eta * incident + (eta * cosi - sqrtf(k)) * n;
}



RayTriangleIntersection Get_closest_intersection(vec3 cameraPosition, vec3 rayDirection, const std::vector<ModelTriangle>& triangles,int depth) {
    RayTriangleIntersection closestIntersection;

    closestIntersection.distanceFromCamera = numeric_limits<float>::infinity();
    for(int i = 0; i < triangles.size(); i++) {
        ModelTriangle triangle = triangles[i];

        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 spVector = cameraPosition - triangle.vertices[0];
        mat3 DEMatrix(-rayDirection, e0, e1);
        vec3 possibleSolution = inverse(DEMatrix) * spVector;

        float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
        float w =1-u-v;
        glm::vec2 textureCoords = w *glm::vec2(triangle.texturePoints[0].x,triangle.texturePoints[0].y)+
                                  u *glm::vec2(triangle.texturePoints[1].x,triangle.texturePoints[1].y)+
                                  v *glm::vec2(triangle.texturePoints[2].x,triangle.texturePoints[2].y);
        if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0 && t > 0) {
            if(t < closestIntersection.distanceFromCamera) {
                vec3 intersectionPoint=cameraPosition+possibleSolution.x*rayDirection;
                closestIntersection.distanceFromCamera = t;
                closestIntersection.triangleIndex = i;
                closestIntersection.u = u;
                closestIntersection.v = v;
                closestIntersection=RayTriangleIntersection(intersectionPoint,closestIntersection.distanceFromCamera,triangle,closestIntersection.triangleIndex,u,v);
                closestIntersection.Texture=textureCoords;
            }
        }
    }
    if(depth>=5){
        closestIntersection.intersectedTriangle.colour = Colour(255,255,255);
        return closestIntersection;
    }
    if(closestIntersection.intersectedTriangle.mirror && depth<5){

        vec3 normal= calculateFaceNormal(closestIntersection.intersectedTriangle);
        glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
        RayTriangleIntersection mirrorIntersection=Get_closest_intersection(closestIntersection.intersectionPoint + reflectedDirection * 0.001f,reflectedDirection,triangles,depth+1);
        if(!isinf(mirrorIntersection.distanceFromCamera)){
            closestIntersection=mirrorIntersection;
        }

    }else if(closestIntersection.intersectedTriangle.glass && depth<5){

        vec3 normal= calculateFaceNormal(closestIntersection.intersectedTriangle);
        float refractiveIndex = 1.5;

        glm::vec3 reflectedDirection = refract(rayDirection, normal, refractiveIndex);

        RayTriangleIntersection glassIntersection=Get_closest_intersection(closestIntersection.intersectionPoint + reflectedDirection * 0.001f,reflectedDirection,triangles,depth+1);
        if(!isinf(glassIntersection.distanceFromCamera)){
            closestIntersection=glassIntersection;
        }

    }
    else if(closestIntersection.intersectedTriangle.metal && depth<5){
        glm::vec3 normal = calculateFaceNormal(closestIntersection.intersectedTriangle);
        glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
        RayTriangleIntersection metalIntersection = Get_closest_intersection(closestIntersection.intersectionPoint + reflectedDirection * 0.001f, reflectedDirection, triangles, depth+1);
        if (!isinf(metalIntersection.distanceFromCamera)) {
            closestIntersection = metalIntersection;
        }
    }

    return closestIntersection;
}


int IsShadow(RayTriangleIntersection intersection, vector<vec3>& lightSources, const vector<ModelTriangle>& triangles) {
    int unblockedSamples = 0;
    for (const auto& lightPosition : lightSources) {
        vec3 shadow_ray = normalize(lightPosition - intersection.intersectionPoint);
        bool blocked = false;
        for (int i = 0; i < triangles.size(); i++) {
            ModelTriangle triangle = triangles[i];

            vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
            vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
            vec3 spVector = intersection.intersectionPoint - triangle.vertices[0];
            mat3 DEMatrix(-shadow_ray, e0, e1);
            vec3 possibleSolution = inverse(DEMatrix) * spVector;

            float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
            if (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0 && t > 0.01) {
                if (t < length(lightPosition - intersection.intersectionPoint) && i != intersection.triangleIndex) {
                    blocked = true;
                    break;
                }
            }
        }
        if (!blocked) {
            unblockedSamples++;
        }
    }
    return unblockedSamples;
}







float getProximityLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint) {
    float distance = glm::length(lightPosition - intersectionPoint);
    return std::min(1.0f, std::max(0.0f, 20.0f / static_cast<float>(4 * M_PI * distance * distance)));
}

float getIncidenceLighting(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint, const glm::vec3& normal) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    return std::max(0.0f, glm::dot(glm::normalize(normal), lightDirection));
    //return glm::clamp(glm::dot(glm::normalize(normal), lightDirection), 0.0f, 1.0f);
}

float getSpecularIntensity(const glm::vec3& lightPosition, const glm::vec3& intersectionPoint,  const glm::vec3& normal, float specularExponent,const glm::vec3 cameraPosition,mat3 Camera_Orientation) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - intersectionPoint);
    glm::vec3 viewDirection = glm::normalize(cameraPosition-intersectionPoint);
    //viewDirection=glm::normalize(Camera_Orientation*viewDirection);
    glm::vec3 reflectedRay = glm::reflect(-lightDirection, normal);
    //glm::vec3 reflectedRay = (-lightDirection) - (2.0f * normal) * (glm::dot(glm::normalize(-lightDirection), glm::normalize(normal)));
    // float reflectionAngle = glm::dot(glm::normalize(-lightDirection), glm::normalize(reflectedRay));
    float reflectionAngle = glm::dot(glm::normalize(viewDirection), glm::normalize(reflectedRay));
    reflectionAngle=glm::max(0.0f,reflectionAngle);
    float SpecularLighting = pow(reflectionAngle, specularExponent);
    return SpecularLighting;
}






float getBrightness( RayTriangleIntersection &intersection, glm::vec3 &normal, glm::vec3 cameraPosition, float specularExponent, const vector<ModelTriangle>& triangles, mat3 Camera_Orientation, bool diffuse, bool specular){
    float totalDiffuseIntensity = 0.0f;
    float maxSpecularIntensity = 0.0f;
    std::vector<glm::vec3> lightSource = lightSources();

    for (const auto& lightPosition : lightSource) {
        float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
        float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, normal);
        totalDiffuseIntensity += ProximityLighting * IncidenceLighting;
        if (specular) {
            float specularLighting = getSpecularIntensity(lightPosition, intersection.intersectionPoint, normal, specularExponent, cameraPosition, Camera_Orientation);
            maxSpecularIntensity = std::max(maxSpecularIntensity, specularLighting);
        }
    }

    float ambientIntensity = 0.1f;
    float totalBrightness;
    if(diffuse){
        totalBrightness = ambientIntensity + totalDiffuseIntensity / lightSource.size();
    }

    if (specular) {
        totalBrightness = ambientIntensity + totalDiffuseIntensity / lightSource.size()+maxSpecularIntensity;
    }

    return glm::clamp(totalBrightness, 0.2f, 1.0f);
}


float getVertexBrightness(const RayTriangleIntersection &intersection, vec3& lightPosition, const glm::vec3 &normal, const glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation,const vector<ModelTriangle>& triangles) {
    float ProximityLighting = getProximityLighting(lightPosition, intersection.intersectionPoint);
    //incidence
    float IncidenceLighting = getIncidenceLighting(lightPosition, intersection.intersectionPoint, normal);

    float specularLighting = getSpecularIntensity(lightPosition, intersection.intersectionPoint, normal,specularExponent, cameraPosition, Camera_Orientation);
    //get brightness
    float ambientIntensity = 0.1f;

    float brightness = ambientIntensity + ProximityLighting * IncidenceLighting + specularLighting;

    brightness = glm::clamp(brightness, 0.0f, 1.0f);

    return brightness;
}


vector<ModelTriangle> getVertexNormals(vector<ModelTriangle> &triangles) {
    for ( ModelTriangle& triangle : triangles) {
        glm::vec3 normal = calculateFaceNormal(triangle);
        triangle.normal = glm::normalize(normal);
    }

    for (ModelTriangle& triangle : triangles) {
        for (int i = 0; i < 3; ++i) {
            glm::vec3 vertexNormal = glm::vec3(0.0f);
            int count = 0;
            // Find triangles that share the same vertices
            for (ModelTriangle& sharedTriangle : triangles) {
                glm::vec3 vertex = triangle.vertices[i];
                if (vertex == sharedTriangle.vertices[0]||vertex == sharedTriangle.vertices[1]||vertex == sharedTriangle.vertices[2]) {
                    vertexNormal += sharedTriangle.normal;
                    count++;
                }
            }
            if (count > 0) {
                glm::vec3 averageNormal = glm::normalize(vertexNormal / static_cast<float>(count));
                triangle.vertexNormals[i] = averageNormal;

            }
        }
    }
    return triangles;
}


glm::vec3 getInterpolatedNormal(const ModelTriangle &triangle, float u, float v) {
    float w = 1.0f - u - v;
    glm::vec3 interpolatedNormal = w * triangle.vertexNormals[0] + u * triangle.vertexNormals[1] + v * triangle.vertexNormals[2];
    return glm::normalize(interpolatedNormal);
}


Colour gouraudShading(const RayTriangleIntersection &intersection,  glm::vec3 &lightPosition, const glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation,const vector<ModelTriangle>& triangles) {
    const ModelTriangle &triangle = intersection.intersectedTriangle;
    std::vector<float> vertexLightIntensities(3);
    for (int i = 0; i < 3; ++i) {
        glm::vec3 vertex = triangle.vertices[i];
        glm::vec3 vertexNormal = triangle.vertexNormals[i];

        float vertexBrightness = getVertexBrightness(intersection, lightPosition, vertexNormal, cameraPosition, specularExponent,Camera_Orientation,triangles);
        vertexLightIntensities[i] = vertexBrightness;
    }

    float u = intersection.u;
    float v = intersection.v;
    float w = 1.0f - u - v;
    float interpolatedLightIntensity = w * vertexLightIntensities[0] + u * vertexLightIntensities[1] + v * vertexLightIntensities[2];

    interpolatedLightIntensity=glm::clamp(interpolatedLightIntensity,0.0f,1.0f);
    Colour colour = intersection.intersectedTriangle.colour;
    colour.red = static_cast<uint8_t>(colour.red * interpolatedLightIntensity);
    colour.green = static_cast<uint8_t>(colour.green * interpolatedLightIntensity);
    colour.blue = static_cast<uint8_t>(colour.blue * interpolatedLightIntensity);

    return colour;
}


Colour phongShading(RayTriangleIntersection &intersection, glm::vec3 &lightPosition, glm::vec3 &cameraPosition, float specularExponent,mat3 Camera_Orientation) {
    glm::vec3 interpolatedNormal = getInterpolatedNormal(intersection.intersectedTriangle, intersection.u, intersection.v);

    float ProximityLighting= getProximityLighting(lightPosition,intersection.intersectionPoint);
    float IncidenceLighting   = getIncidenceLighting(lightPosition,intersection.intersectionPoint,interpolatedNormal);
    float specularLighting = getSpecularIntensity(lightPosition,intersection.intersectionPoint,interpolatedNormal,specularExponent,cameraPosition,Camera_Orientation);

    Colour colour = intersection.intersectedTriangle.colour;
    colour.red = static_cast<uint8_t>(std::min(colour.red * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
    colour.green = static_cast<uint8_t>(std::min(colour.green * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));
    colour.blue = static_cast<uint8_t>(std::min(colour.blue * ProximityLighting*IncidenceLighting + specularLighting * 255, 255.0f));

    return colour;
}



glm::vec3 getRayDirection(float x,float y,glm::vec3 camera_position,mat3 Camera_Orientation){
    float scalingFactor = 150.0f;
    float focal_length = 2.0f;

    double x1 = (x- camera_position.x - (WIDTH / 2)) / (scalingFactor * focal_length);
    double y1 = -(y-camera_position.y - (HEIGHT / 2)) / (scalingFactor * focal_length);
    double z1 = -focal_length ;

    glm::vec3 ray_direction = glm::vec3(x1 , y1 , z1);
    glm::vec3 direction= normalize(Camera_Orientation*ray_direction);
    return direction;
}



glm::vec3 randomInHemisphere(const glm::vec3 &normal) {
    glm::vec3 inUnitSphere;
    do {
        inUnitSphere = 2.0f * glm::vec3(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX) - glm::vec3(1, 1, 1);
    } while (glm::length(inUnitSphere) >= 1.0f);

    if (glm::dot(inUnitSphere, normal) < 0.0f) {
        inUnitSphere = -inUnitSphere;
    }

    return inUnitSphere;
}


Colour CalculateColor(const vector<ModelTriangle>& triangles,RayTriangleIntersection &intersection, glm::vec3 &cameraPosition, float specularExponent, mat3 &Camera_Orientation,  glm::vec3 rayDirection, int depth=0) {
    if (depth >= 5) return Colour(0,0,0);
    Colour colour;
    vector<vec3> lightSource=lightSources();
    vector<ModelTriangle> sphere1;
    vector<ModelTriangle> sphere2;
    vector<ModelTriangle> sphere3;

    vector<ModelTriangle> normalTriangles;



    glm::vec3 normal = calculateNormal(intersection.intersectedTriangle, intersection);
    int unblockedSamples = IsShadow(intersection, lightSource, triangles);
    float shadowIntensity = static_cast<float>(unblockedSamples) / static_cast<float>(lightSource.size());

    for(auto& triangle:triangles){
        if(triangle.material=="sphere1"){
            sphere1.push_back(triangle);
        }else if(triangle.material=="sphere2"){
            sphere2.push_back(triangle);
        }else if(triangle.material=="sphere3"){
            sphere3.push_back(triangle);
        }else{
            normalTriangles.push_back(triangle);
        }
    }

    TextureMap textureFile;




    if(intersection.intersectedTriangle.material=="sphere3"){
        colour = gouraudShading(intersection, lightPosition, cameraPosition, specularExponent, Camera_Orientation, triangles);
        //cout<<"g:"<<intersection.intersectedTriangle.material<<endl;
        //cout<<"g:"<<colour<<endl;

    }
    else if(intersection.intersectedTriangle.material=="sphere1"){
        colour = phongShading(intersection, lightPosition, cameraPosition, specularExponent, Camera_Orientation);
        //cout<<"p:"<<intersection.intersectedTriangle.material<<endl;

    }
    else if(intersection.intersectedTriangle.material=="sphere2"){
        //specular
        float brightness = getBrightness(intersection, normal, cameraPosition, specularExponent,sphere3,Camera_Orientation,false,true);
        shadowIntensity=clamp(shadowIntensity,0.4f,1.0f);
        brightness *= shadowIntensity;
        colour = intersection.intersectedTriangle.colour;
        colour.red *= std::min(brightness, 1.0f);
        colour.green *= std::min(brightness, 1.0f);
        colour.blue *= std::min(brightness, 1.0f);
        //cout<<"g:"<<colour<<endl;

    }
    else  {
        float brightness = getBrightness(intersection, normal, cameraPosition, specularExponent, normalTriangles,Camera_Orientation, true, false);
        shadowIntensity = clamp(shadowIntensity, 0.4f, 1.0f);
        brightness *= shadowIntensity;
        colour = intersection.intersectedTriangle.colour;
        colour.red *= std::min(brightness, 1.0f);
        colour.green *= std::min(brightness, 1.0f);
        colour.blue *= std::min(brightness, 1.0f);
    }

    if (intersection.intersectedTriangle.mirror&&depth<5) {
        glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
        RayTriangleIntersection reflectedIntersection = Get_closest_intersection(intersection.intersectionPoint + reflectedDirection * 0.001f, reflectedDirection, triangles,depth+1);
        if (!isinf(intersection.distanceFromCamera)) {
            Colour reflectedColour = CalculateColor(triangles, reflectedIntersection, cameraPosition,
                                                    specularExponent, Camera_Orientation, rayDirection,depth+1);
            colour = reflectedColour;


        }

    }else if(intersection.intersectedTriangle.glass&&depth<5){

        float cosi = glm::clamp(glm::dot(rayDirection, normal), -1.0f, 1.0f);
        float etai = 1, etat = 1.5f;
        if (cosi > 0) { std::swap(etai, etat); }
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        float reflectivity;
        if (sint >= 1) {
            reflectivity = 1;
        } else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            reflectivity = (Rs * Rs + Rp * Rp) / 2;
        }


        glm::vec3 reflectedDirection = glm::normalize(rayDirection) - 2 * glm::dot(rayDirection, normal) * normal;
        glm::vec3 refractedRay = refract(rayDirection, normal, etai / etat);


        RayTriangleIntersection reflectedIntersection = Get_closest_intersection(intersection.intersectionPoint + reflectedDirection * 0.001f, reflectedDirection, triangles, depth+1);
        RayTriangleIntersection refractionIntersection = Get_closest_intersection(intersection.intersectionPoint + refractedRay * 0.001f, refractedRay, triangles, depth+1);

        Colour reflectedColour = CalculateColor(triangles, reflectedIntersection, cameraPosition, specularExponent, Camera_Orientation, reflectedDirection, depth+1);
        Colour refractionColour = CalculateColor(triangles, refractionIntersection, cameraPosition, specularExponent, Camera_Orientation, refractedRay, depth+1);


        colour.red = reflectedColour.red * reflectivity + refractionColour.red * (1 - reflectivity);
        colour.green = reflectedColour.green * reflectivity + refractionColour.green * (1 - reflectivity);
        colour.blue = reflectedColour.blue * reflectivity + refractionColour.blue * (1 - reflectivity);

    }else if(intersection.intersectedTriangle.texture&&depth<5){
        textureFile=TextureFile1;
        float u = intersection.u;
        float v = intersection.v;
        float w = 1.0f - u - v;
        float textureX = (w * intersection.intersectedTriangle.texturePoints[0].x + u * intersection.intersectedTriangle.texturePoints[1].x + v * intersection.intersectedTriangle.texturePoints[2].x) * textureFile.width;
        float textureY = (w * intersection.intersectedTriangle.texturePoints[0].y + u * intersection.intersectedTriangle.texturePoints[1].y + v * intersection.intersectedTriangle.texturePoints[2].y) * textureFile.height;


        textureX = std::fmod(textureX, textureFile.width);
        textureY = std::fmod(textureY, textureFile.height);


        uint32_t colour_val = textureFile.pixels[static_cast<int>(textureY) * textureFile.width + static_cast<int>(textureX)];
        int red = (colour_val >> 16) & 0xFF;
        int green = (colour_val >> 8) & 0xFF;
        int blue = colour_val & 0xFF;


        colour.red = red;
        colour.green = green;
        colour.blue = blue;


    }
    else if(intersection.intersectedTriangle.metal&&depth<5){

        glm::vec3 normal = calculateFaceNormal(intersection.intersectedTriangle);
        glm::vec3 reflectedDirection = glm::normalize(rayDirection - 2 * glm::dot(rayDirection, normal) * normal);
        float roughnessFactor = 0.05;
        glm::vec3 randomVector = randomInHemisphere(normal);
        reflectedDirection += roughnessFactor * randomVector;
        reflectedDirection = glm::normalize(reflectedDirection);
        RayTriangleIntersection metalIntersection = Get_closest_intersection(intersection.intersectionPoint + reflectedDirection * 0.001f, reflectedDirection, triangles, depth + 1);
        Colour reflectedColour;
        if (!isinf(metalIntersection.distanceFromCamera)) {
            reflectedColour = CalculateColor(triangles, metalIntersection, cameraPosition, specularExponent, Camera_Orientation, reflectedDirection, depth + 1);

        }
        //colour = reflectedColour;
        float reflectivity = 0.5;
        Colour goldBaseColour(255, 215, 0);

        colour.red = glm::mix(goldBaseColour.red, reflectedColour.red, reflectivity);
        colour.green = glm::mix(goldBaseColour.green, reflectedColour.green, reflectivity);
        colour.blue = glm::mix(goldBaseColour.blue, reflectedColour.blue, reflectivity);
        //cout<<"metal"<<colour<<endl;
    }
    else if(intersection.intersectedTriangle.texture){
        TextureMap& textureFile=TextureFile1;

        // cout<<"1"<<endl;

        float u = intersection.u;
        float v = intersection.v;
        float w = 1.0f - u - v;
        float textureX = (w * intersection.intersectedTriangle.texturePoints[0].x + u * intersection.intersectedTriangle.texturePoints[1].x + v * intersection.intersectedTriangle.texturePoints[2].x) * textureFile.width;
        float textureY = (w * intersection.intersectedTriangle.texturePoints[0].y + u * intersection.intersectedTriangle.texturePoints[1].y + v * intersection.intersectedTriangle.texturePoints[2].y) * textureFile.height;


        textureX = std::fmod(textureX, textureFile.width);
        textureY = std::fmod(textureY, textureFile.height);


        uint32_t colour_val = textureFile.pixels[static_cast<int>(textureY) * textureFile.width + static_cast<int>(textureX)];
        int red = (colour_val >> 16) & 0xFF;
        int green = (colour_val >> 8) & 0xFF;
        int blue = colour_val & 0xFF;


        colour.red = red;
        colour.green = green;
        colour.blue = blue;


    }



    return colour;
}




void DrawRay(const std::vector<ModelTriangle>& triangles, DrawingWindow &window, glm::vec3 camera_position,mat3 Camera_Orientation,float specularExponent) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            glm::vec3 rayDirection= getRayDirection(x,y,camera_position,Camera_Orientation);
            RayTriangleIntersection intersection = Get_closest_intersection(camera_position, rayDirection, triangles,1);
            if (!isinf(intersection.distanceFromCamera)) {
                Colour colour = CalculateColor(triangles,intersection, camera_position, specularExponent, Camera_Orientation, rayDirection);

                uint32_t color_value =
                        (255 << 24) + (static_cast<int>(colour.red) << 16) + (static_cast<int>(colour.green) << 8) +
                        static_cast<int>(colour.blue);
                window.setPixelColour(x, y, color_value);
            }

        }
    }
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
    std::vector<ModelTriangle> modelTriangles1=readOBJ("cornell-box.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles2=readOBJ("textured-cornell-box.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles3=readOBJ("sphere.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles4=readOBJ("loveBox3.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles5=readOBJ("loveBox.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles6=readOBJ("cornell-box2.obj", 0.35);
    std::vector<ModelTriangle> modelTriangles7=readOBJ("cornell-box3.obj", 0.35);



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
        modelTriangles6= getVertexNormals(modelTriangles6);

        //normal wireframe
        //RenderScene(window, modelTriangles1, cameraPosition, Camera_Orientation, depthBuffer, true);

        //normal render
        RenderScene(window, modelTriangles1, cameraPosition, Camera_Orientation, depthBuffer, false);

        //render texture
        //RenderScene(window, modelTriangles2, cameraPosition, Camera_Orientation, depthBuffer, false);

        //raytrace texture
        //DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);

        // Metal
        //DrawRay(modelTriangles7, window, cameraPosition, Camera_Orientation,  256.0f);

        //Glass
        //DrawRay(modelTriangles6, window, cameraPosition, Camera_Orientation,  256.0f);

        //normal raytrace
        //DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);

        //love box or using Key 5
        //DrawRay(modelTriangles4, window, cameraPosition, Camera_Orientation, 256.0f);

        //vec3 lightPosition(0,0,0);
        if (model == WireFrame) {
            RenderScene(window, modelTriangles1, cameraPosition, Camera_Orientation, depthBuffer, true);
        }
        else if (model == Render) {
            RenderScene(window, modelTriangles2, cameraPosition, Camera_Orientation, depthBuffer, false);
        }
        else if (model == Diffuse) {
            DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);
        }
        else if (model == Texture) {
            DrawRay(modelTriangles5, window, cameraPosition, Camera_Orientation,  256.0f);
        }
        else if (model == Specular) {
            DrawRay(modelTriangles4, window, cameraPosition, Camera_Orientation, 256.0f);
        }
        else if (model == Soft) {
            DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
        }
        else if (model == Mirror) {
            DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);
        }

//
//        for (int frameCount = 0; frameCount < 10; frameCount++) {
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


//        for (int frameCount = 0; frameCount < 10; frameCount++) {
//            window.clearPixels();
//            int n_zero = 5;
//            //DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);
//            //DrawRay(modelTriangles4, window, cameraPosition, Camera_Orientation, 256.0f);
//            DrawRay(modelTriangles6, window, cameraPosition, Camera_Orientation,  256.0f);
//            //DrawRay(modelTriangles2, window, cameraPosition, Camera_Orientation,  256.0f);
//            //DrawRay(modelTriangles1, window, cameraPosition, Camera_Orientation,  256.0f);
//
//            //RenderScene(window, modelTriangles4, cameraPosition, Camera_Orientation, depthBuffer, true);
//            string name = string(n_zero - to_string(frameCount).length(), '0') + to_string(frameCount);
//            window.savePPM("../output1/"+name+".ppm");
//            cout << "saved " << frameCount << endl;
////
////           if(frameCount < 43) {
////               cameraPosition.z-=0.06f;
////           }
////        else if(frameCount < 49) {
////               cameraPosition.x+=0.06f;
////
////           }
////        else if(frameCount < 56) {
////               cameraPosition.y-=0.06f;
//
////           }else if(frameCount < 41) {
////               cameraPosition.x-=0.06f;
////
////           }else if(frameCount < 48) {
////               cameraPosition.y-=0.06f;
////
//           //}
//
//
//        }
//
//animationCompleted = true;
//





        window.renderFrame();



    }


}