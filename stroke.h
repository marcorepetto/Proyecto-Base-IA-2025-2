#ifndef STROKE_H
#define STROKE_H

#include <vector>
#include <string>
#include <cstdint>

// ================= Canvas =================
struct Canvas {
    int width, height;
    std::vector<float> rgb; // tamaño = width * height * 3

    Canvas(int w, int h);
    void clear(float r = 255, float g = 255, float b = 255);
    void setPixel(int x, int y, float r, float g, float b);
    std::vector<uint8_t> getRGBData() const;
};

// Imagen en escala de grises (para brushes)
struct ImageGray {
    int width = 0, height = 0;
    std::vector<float> gray; // tamaño = width * height
};

bool loadImageGray(const std::string& filename, ImageGray& out);
bool savePNG(const Canvas& C, const std::string& filename);

// ================= Stroke =================
class Stroke {
public:
    float x_rel = 0.5f;
    float y_rel = 0.5f;
    float size_rel = 0.2f;
    float rotation_deg = 0.0f;
    int   type = 0;
    float r = 0, g = 0, b = 0;
    std::vector<float> strokeAlphas;

    Stroke() = default;
    Stroke(float xr, float yr, float sr, float rot, int t,
           float rr, float gg, float bb, int canvas_width = 48, int canvas_height = 64);

    void draw(Canvas C);
    void randomize();
};


// Pinta el lienzo con todos los strokes (en el orden recibido)
bool loadImageRGB_asCanvas(const std::string& filename, Canvas& out);

extern std::vector<ImageGray> gBrushes;

#endif
