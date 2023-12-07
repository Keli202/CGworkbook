#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Utils.h"

class TextureMap {
public:
	size_t width;
	size_t height;
	std::vector<uint32_t> pixels;


    uint32_t getPixel(int x, int y) const {
        if (x < 0 || x >= width || y < 0 || y >= height) {
            return 0;
        }

        return pixels[y * width + x];
    }
	TextureMap();
	TextureMap(const std::string &filename);
	friend std::ostream &operator<<(std::ostream &os, const TextureMap &point);
};
