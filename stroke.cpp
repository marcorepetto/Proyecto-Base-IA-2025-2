#include "stroke.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <chrono>

#include <Eigen/Dense>

// stb_image / stb_image_write
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

std::vector<ImageGray> gBrushes;  // definición del vector global


// Canvas

Canvas::Canvas(int w, int h) : width(w), height(h), rgb(w*h*3, 0.0f) {} // arranca blanco

void Canvas::clear(float r, float g, float b) {
    const int N = width * height;
    for (int i = 0; i < N; i++) {
        rgb[i*3 + 0] = r;
        rgb[i*3 + 1] = g;
        rgb[i*3 + 2] = b;
    }
}

void Canvas::setPixel(int x, int y, float r, float g, float b) {
    if (x < 0 || y < 0 || x >= width || y >= height) return;
    const int idx = (y * width + x) * 3;
    rgb[idx + 0] = r;
    rgb[idx + 1] = g;
    rgb[idx + 2] = b;
}

std::vector<uint8_t> Canvas::getRGBData() const {
    std::vector<uint8_t> data(width * height * 3);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = static_cast<uint8_t>(std::clamp(rgb[i]*255.0f, 0.0f, 255.0f));
    }
    return data;
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
    out.gray.assign(data, data + (w*h));
    stbi_image_free(data);
    return true;
}

// Helpers para guardar la imagen final

bool savePNG(const Canvas& C, const std::string& filename) {
    const int stride = C.width * 3;
    if (!stbi_write_png(filename.c_str(), C.width, C.height, 3, C.getRGBData().data(), stride)) {
        std::cerr << "stb_image_write: no pude guardar " << filename << "\n";
        return false;
    }
    return true;
}

// Helpers para comparar lienzos

bool loadImageRGB_asCanvas(const std::string& filename, Canvas& out) {
    int w, h, n;
    unsigned char* data = stbi_load(filename.c_str(), &w, &h, &n, 3);
    if (!data) {
        std::cerr << "stb_image: no pude cargar " << filename << "\n";
        return false;
    }
    out = Canvas(w, h);
    for (int i = 0; i < w * h * 3; ++i) {
        float pixel = data[i] / 255.0f;
        out.rgb[i] = pixel;
    }
    stbi_image_free(data);
    return true;
}

// Stroke
Stroke::Stroke(float xr, float yr, float sr, float rot, int t,
               float rr, float gg, float bb, int canvas_width, int canvas_height)
    : x_rel(xr), y_rel(yr), size_rel(sr), rotation_deg(rot),
      type(t), r(rr), g(gg), b(bb) {
        strokeAlphas.resize(canvas_width * canvas_height, 0.0f);
      }

template <typename T>
static inline T clampT(T v, T lo, T hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

void Stroke::draw(Canvas C) {
    // --- Validaciones básicas ---
    if (gBrushes.empty()) return;
    if (type < 0 || type >= (int)gBrushes.size()) return;

    const ImageGray& brush = gBrushes[type];
    if (brush.width == 0 || brush.height == 0) return;

    // Limpiar alphas previos
    std::fill(strokeAlphas.begin(), strokeAlphas.end(), 0.0f);

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
                return brush.gray[y * bw + x] / 255.0f;  // [0,1]
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

            // 6) GUARDAR ALPHAS
            strokeAlphas[y_dst * C.width + x_dst] = a;
        }
    }
}

void Stroke::randomize() {
    // Random engine
    static std::random_device rd;
    static std::mt19937 gen(rd());

    // Random distributions
    std::uniform_real_distribution<float> dist_pos(0.0f, 1.0f);
    std::uniform_real_distribution<float> dist_size(0.15f, 0.7f);
    std::uniform_real_distribution<float> dist_rot(0.0f, 360.0f);
    std::uniform_int_distribution<int>    dist_type(0, (int)gBrushes.size() - 1);
    std::uniform_real_distribution<float> dist_color(0.0f, 1.0f);

    // Randomize attributes
    x_rel = dist_pos(gen);
    y_rel = dist_pos(gen);
    size_rel = dist_size(gen);
    rotation_deg = dist_rot(gen);
    type = dist_type(gen);

    // std::cout << "Randomized Stroke: "
    //      << "x_rel=" << x_rel << ", y_rel=" << y_rel
    //      << ", size_rel=" << size_rel << ", rotation_deg=" << rotation_deg
    //      << ", type=" << type << "\n";
}