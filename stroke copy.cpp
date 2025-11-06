#include "stroke.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

// stb_image / stb_image_write
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

std::vector<ImageGray> gBrushes;  // definición del vector global

// Canvas

Canvas::Canvas(int w, int h) : width(w), height(h), rgb(w*h*3, 255) {} // arranca blanco

void Canvas::clear(uint8_t r, uint8_t g, uint8_t b) {
    const int N = width * height;
    for (int i = 0; i < N; i++) {
        rgb[i*3 + 0] = r;
        rgb[i*3 + 1] = g;
        rgb[i*3 + 2] = b;
    }
}

void Canvas::setPixel(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
    if (x < 0 || y < 0 || x >= width || y >= height) return;
    const int idx = (y * width + x) * 3;
    rgb[idx + 0] = r;
    rgb[idx + 1] = g;
    rgb[idx + 2] = b;
}

// Helper para cargar los brushes

bool loadImageGray(const std::string& filename, ImageGray& out) {
    int w, h, n;
    unsigned char* data = stbi_load(filename.c_str(), &w, &h, &n, 1);
    if (!data) {
        std::cerr << "stb_image: no pude cargar " << filename << "\n";
        return false;
    }
    out.width = w;
    out.height = h;
    out.data.assign(data, data + (w*h));
    stbi_image_free(data);
    return true;
}

// Helpers para guardar la imagen final

bool savePNG(const Canvas& C, const std::string& filename) {
    const int stride = C.width * 3;
    if (!stbi_write_png(filename.c_str(), C.width, C.height, 3,
                        C.rgb.data(), stride)) {
        std::cerr << "stb_image_write: no pude guardar " << filename << "\n";
        return false;
    }
    return true;
}

// Helpers para comparar lienzos

void render(const std::vector<Stroke>& strokes, Canvas& C) {
    C.clear(255, 255, 255); // fondo blanco
    for (const auto& s : strokes) s.draw(C);
}

bool loadImageRGB_asCanvas(const std::string& filename, Canvas& out) {
    int w, h, n;
    unsigned char* data = stbi_load(filename.c_str(), &w, &h, &n, 3);
    if (!data) {
        std::cerr << "stb_image: no pude cargar " << filename << "\n";
        return false;
    }
    out = Canvas(w, h);
    out.rgb.assign(data, data + (w*h*3));
    stbi_image_free(data);
    return true;
}

// Stroke

Stroke::Stroke(float xr, float yr, float sr, float rot, int t,
               uint8_t rr, uint8_t gg, uint8_t bb)
    : x_rel(xr), y_rel(yr), size_rel(sr), rotation_deg(rot),
      type(t), r(rr), g(gg), b(bb) {}

template <typename T>
static inline T clampT(T v, T lo, T hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

void Stroke::draw(Canvas& C) const {
    // --- Validaciones básicas ---
    if (gBrushes.empty()) return;
    if (type < 0 || type >= (int)gBrushes.size()) return;

    const ImageGray& brush = gBrushes[type];
    if (brush.width == 0 || brush.height == 0) return;

    // 1) CONFIGURACIÓN DE TAMAÑO / ESCALA (mantener aspecto del brush)
    //    - 'base' controla el tamaño general
    const int base = std::max(1, int(size_rel * std::min(C.width, C.height)));
    const int bw = brush.width;
    const int bh = brush.height;

    // Escala para que el lado MAYOR del brush quede en 'base'
    const float s    = float(base) / float(std::max(bw, bh));  
    const float invs = (s > 0.0f) ? (1.0f / s) : 0.0f;          

    // Tamaño real que ocupará en el canvas respetando el aspecto
    const int w_pix = std::max(1, int(std::round(bw * s)));
    const int h_pix = std::max(1, int(std::round(bh * s)));
    const int halfW = w_pix / 2;
    const int halfH = h_pix / 2;

    // 2) TRASLACIÓN (MOVER): centro de la pincelada en el canvas
    const int cx = clampT(int(std::round(x_rel * C.width)),  0, C.width  - 1);
    const int cy = clampT(int(std::round(y_rel * C.height)), 0, C.height - 1);

    // 3) ROTACIÓN: precomputar cos/sin del ángulo (en radianes)
    const float PI = 3.14159265358979323846f;
    const float theta = rotation_deg * (PI / 180.0f);
    const float ct = std::cos(theta);
    const float st = std::sin(theta);

    // 4) RASTERIZADO: recorrer el rectángulo destino w_pix x h_pix
    //    Para cada píxel destino, aplicamos la TRANSFORMACIÓN INVERSA:
    //      - des-rotar (R(-θ))
    //        R(-θ) = | cosθ   sinθ |
    //                | -sinθ  cosθ |
    //      - des-escalar (S(1/s))
    //    para obtener la coordenada (xb,yb) en el espacio del brush.
    //    Luego muestreamos con bilineal y mezclamos color.
    for (int dy = -halfH; dy <= halfH; ++dy) {
        for (int dx = -halfW; dx <= halfW; ++dx) {

            // -------- LÍMITES DEL CANVAS --------
            const int x_dst = cx + dx;   // TRASLACIÓN aplicada aquí
            const int y_dst = cy + dy;
            if (x_dst < 0 || y_dst < 0 || x_dst >= C.width || y_dst >= C.height) continue;

            // -------- ROTACIÓN (inversa) --------
            // Vector relativo al centro (dx,dy) -> des-rotado
            const float xr =  dx * ct + dy * st;   // R(-θ) * [dx, dy]
            const float yr = -dx * st + dy * ct;

            // -------- ESCALA (inversa) --------
            // Pasar a coordenadas del brush en píxeles, centrado en (0,0)
            const float xb = xr * invs;
            const float yb = yr * invs;

            // Si cae fuera del brush original, omitimos
            if (std::fabs(xb) > (bw - 1) * 0.5f || std::fabs(yb) > (bh - 1) * 0.5f) continue;

            // -------- COORDS NORMALIZADAS [0,1] PARA MUESTREO --------
            const float tu = (xb / float(bw - 1)) + 0.5f;
            const float tv = (yb / float(bh - 1)) + 0.5f;

            // 5) MUESTREO BILINEAL DEL BRUSH (canal "alpha"/máscara)
            const float fu = tu * (bw - 1);
            const float fv = tv * (bh - 1);
            int   x0 = (int)std::floor(fu);
            int   y0 = (int)std::floor(fv);
            int   x1 = clampT(x0 + 1, 0, bw - 1);
            int   y1 = clampT(y0 + 1, 0, bh - 1);
            const float ax = fu - x0;
            const float ay = fv - y0;

            auto sample = [&](int x, int y) -> float {
                return brush.data[y * bw + x] / 255.0f;  // [0,1]
            };
            const float m00 = sample(x0, y0);
            const float m10 = sample(x1, y0);
            const float m01 = sample(x0, y1);
            const float m11 = sample(x1, y1);

            // 'a' = cobertura/alpha del brush en ese punto destino
            const float a = (1 - ax) * (1 - ay) * m00
                          + (    ax) * (1 - ay) * m10
                          + (1 - ax) * (    ay) * m01
                          + (    ax) * (    ay) * m11;

            if (a <= 0.0f) continue;

            // 6) MEZCLA DE COLOR (SRC OVER)
            //    - fg = color del stroke (r,g,b)
            //    - bg = color actual del canvas
            //    - out = a*fg + (1-a)*bg
            const int idx = (y_dst * C.width + x_dst) * 3;

            const float bgR = C.rgb[idx + 0];
            const float bgG = C.rgb[idx + 1];
            const float bgB = C.rgb[idx + 2];

            const float fgR = float(r);
            const float fgG = float(g);
            const float fgB = float(b);

            const float outR = a * fgR + (1.0f - a) * bgR;
            const float outG = a * fgG + (1.0f - a) * bgG;
            const float outB = a * fgB + (1.0f - a) * bgB;

            C.rgb[idx + 0] = (uint8_t)clampT(int(std::round(outR)), 0, 255);
            C.rgb[idx + 1] = (uint8_t)clampT(int(std::round(outG)), 0, 255);
            C.rgb[idx + 2] = (uint8_t)clampT(int(std::round(outB)), 0, 255);
        }
    }
}


