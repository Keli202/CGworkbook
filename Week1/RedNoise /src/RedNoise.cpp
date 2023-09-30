#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>

#define WIDTH 320
#define HEIGHT 240
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

void draw(DrawingWindow &window) {
	window.clearPixels();
    int Width = static_cast<int>(window.width);
    std::vector<float> gradient = interpolateSingleFloats(0.0f, 255.0f, Width);

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




int main(int argc, char *argv[]) {
    std::vector<float> result;
    result = interpolateSingleFloats(2.2, 8.5, 7);
    for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
    std::cout << std::endl;
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

}