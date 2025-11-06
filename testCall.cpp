#include "stroke.h"
#include <iostream>
#include <random>

// RNG global
static std::mt19937 rng(std::random_device{}());
int randInt(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
}

int main() {
    // 1) Cargar brushes desde carpeta
    {
        ImageGray b0, b1, b2, b3;
        if (!loadImageGray("brushes/1.jpg", b0)) return 1;
        if (!loadImageGray("brushes/2.jpg", b1)) return 1;
        if (!loadImageGray("brushes/3.jpg", b2)) return 1;
        if (!loadImageGray("brushes/4.jpg", b3)) return 1;
        gBrushes.push_back(std::move(b0));
        gBrushes.push_back(std::move(b1));
        gBrushes.push_back(std::move(b2));
        gBrushes.push_back(std::move(b3));
    }

    // 2) Canvas, cambiar tamaño de acuerdo a la instancia
    Canvas C(48, 64);

    // 3) Vector de strokes (orden del vector = orden de pintado)
    std::vector<Stroke> strokes;

    // Stroke s (x_rel, y_rel, size_rel, rotation_deg, type, R, G, B)
    Stroke s1(0.50f, 0.50f, 1.35f, 45.0f, 0,randInt(0,255), randInt(0,255), randInt(0,255));
    strokes.push_back(s1);

    Stroke s2(0.70f, 0.40f, 0.95f, 30.0f, 3, randInt(0,255), randInt(0,255), randInt(0,255));
    strokes.push_back(s2);

    Stroke s3(0.20f, 0.20f, 0.5f, 90.0f, 2, randInt(0,255), randInt(0,255), randInt(0,255));
    strokes.push_back(s3);

    // 4) Renderiza (pinta) el vector de strokes en el canva C
    render(strokes, C);

    // 5) Guardar
    if (!savePNG(C, "output.png")) {
        std::cerr << "Error guardando output.png\n";
        return 1;
    }
    std::cout << "OK: guardado output.png\n";
    return 0;
}
