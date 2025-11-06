#include "stroke.h"
#include "model.h"
#include <iostream>
#include <random>


int main(int argc, char** argv) {
    // Leer parámetro 
    //      n: cantidad de brushes a generar
    //      filename: imagen objetivo
    //      n_sample_greedy: cantidad de muestras aleatorias por brush
    //      pixel_threshold: umbral de borde
    //      p: peso de borde respecto a área cubierta
    //      
    int n = std::stoi(argv[1]);
    std::string target_filename = argv[2];
    int n_sample_greedy = std::stoi(argv[3]);
    float pixel_threshold = std::stof(argv[4]);
    float p = std::stof(argv[5]);

    // 1) Cargar brushes desde carpeta
    std::cout << "Cargando brushes...\n";
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

    std::cout << "Brushes cargados: " << gBrushes.size() << "\n";
    // 2) Crear modelo
    std::cout << "Creando modelo con " << n << " strokes...\n";
    Model model = Model(n, target_filename, n_sample_greedy, pixel_threshold, p);

    // 3) Inicializar
    std::cout << "Inicializando...\n";
    model.initialGuess();
    std::cout << "Inicialización completa.\n";
    model.render();

    std::cout << "Calculando pérdida inicial...\n";
    model.current_loss = model.computeLoss();

    std::cout << "Pérdida inicial: " << model.current_loss << "\n";
    std::cout << "Pérdida máxima: " << model.max_loss << "\n";
    std::cout << "Calidad inicial: " << 1.0 - (model.current_loss / model.max_loss) << "\n";

    // 5) Guardar

    // Compose target filename: "n"_"n_sample_greedy"_"pixel_threshold"_"p"_"target_filename"
    std::string output_filename = 
        std::to_string(n) + "_" + 
        std::to_string(n_sample_greedy) + "_" + 
        std::to_string(int(pixel_threshold * 100)) + "_" +
        std::to_string(int(p * 100)) + "_" +
        target_filename;

    if (!savePNG(model.getCurrentCanvas(), "output/" + output_filename)) {
        std::cerr << "Error guardando output/" << output_filename << "\n";
        return 1;
    }
    std::cout << "OK: guardado output/" << output_filename << "\n";
    return 0;
}
