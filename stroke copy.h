#ifndef STROKE_H
#define STROKE_H

#include <vector>
#include <string>
#include <cstdint>

// ================= Canvas =================
struct Canvas {
    int width, height;
    std::vector<uint8_t> rgb; // tamaño = width * height * 3

    Canvas(int w, int h);
    void clear(uint8_t r = 255, uint8_t g = 255, uint8_t b = 255);
    void setPixel(int x, int y, uint8_t r, uint8_t g, uint8_t b);
};

// Imagen en escala de grises (para brushes)
struct ImageGray {
    int width = 0, height = 0;
    std::vector<uint8_t> data;
};

bool loadImageGray(const std::string& filename, ImageGray& out);
bool savePNG(const Canvas& C, const std::string& filename);

// ================= Stroke =================
struct Stroke {
    float x_rel = 0.5f;
    float y_rel = 0.5f;
    float size_rel = 0.2f;
    float rotation_deg = 0.0f;
    int   type = 0;
    uint8_t r = 0, g = 0, b = 0;

    Stroke() = default;
    Stroke(float xr, float yr, float sr, float rot, int t,
           uint8_t rr, uint8_t gg, uint8_t bb);

    void draw(Canvas& C) const;
};

// Pinta el lienzo con todos los strokes (en el orden recibido)
void render(const std::vector<Stroke>& strokes, Canvas& C);

bool loadImageRGB_asCanvas(const std::string& filename, Canvas& out);

// Brushes globales (máscaras en gris). Se cargan en testCall.cpp
extern std::vector<ImageGray> gBrushes;

#endif
