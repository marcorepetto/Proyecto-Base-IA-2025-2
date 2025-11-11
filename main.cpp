#include "stroke.h"
#include "model.h"
#include <iostream>
#include <random>
#include <fstream>
#include <chrono>
#include <cstdlib>

/**
 * @brief Función principal del programa.
 * @param filename Nombre del archivo de la imagen objetivo.
 * @param n Número de strokes a utilizar.
 * @param n_sample_greedy Número de muestras aleatorias por stroke para la inicialización
 * @param p Peso del borde respecto al área cubierta en la inicialización.
 */
int main(int argc, char** argv) {
    std::string target_filename = argv[1];
    int n = std::stoi(argv[2]);
    int n_sample_greedy = std::stoi(argv[3]);
    float p = std::stof(argv[4]);

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
    Model model = Model(n, target_filename, n_sample_greedy, p);

    // 3) Inicialización voraz (timeada)
    auto start = std::chrono::high_resolution_clock::now();
    model.initialGuess();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Inicialización voraz completa. Hecho en " << duration / 1000000.0 << " s.\n";
    // 4) Guardar resultado de inicialización en el archivo output/logs/{model.log_file_name}.txt
    {
        // Crear carpeta de logs
        std::string mkdir_cmd = "mkdir -p output/logs";
        system(mkdir_cmd.c_str());

        std::string log_path = "output/logs/" + model.log_file_name + ".txt";
        std::ofstream logf(log_path);
        if (!logf) {
            std::cerr << "Error creando log '" << log_path << "'\n";
        } else {
            logf << "Model initialization log\n";
            logf << "========================\n";
            logf << "Log file name: " << model.log_file_name << "\n";
            logf << "Target filename: " << target_filename << "\n";
            logf << "Target name (no ext): " << model.target_name << "\n";
            logf << "Number of strokes: " << n << "\n";
            logf << "n_sample_greedy: " << n_sample_greedy << "\n";
            logf << "p: " << p << "\n";
            logf << "Initial guess time (s): " << duration / 1000000.0 << "\n";
            logf << "Initial loss: " << model.current_loss << "\n";
            logf.close();
            std::cout << "OK: guardado log " << log_path << "\n";
        }
    }
    
    // Crear carpeta si no existe
    std::string command = "mkdir -p output/images";
    system(command.c_str());

    // Add model start time to output filename in format targetname_initialGuess_YYYY-MM-DD_HH:MM:SS.png
    std::string output_filename = model.log_file_name + " initialGuess.png";

    if (!savePNG(model.getCurrentCanvas(), "output/images/" + output_filename)) {
        std::cerr << "Error guardando output/" << output_filename << "\n";
        return 1;
    }
    std::cout << "OK: guardado output/" << output_filename << "\n";
    
    {
        std::string log_path = "output/logs/" + model.log_file_name + ".csv";
        std::ofstream logf(log_path);
        if (!logf) {
            std::cerr << "Error creando log '" << log_path << "'\n";
        } else {
            logf << "Iteration" << "," 
                << "Stagnation counter" << "," 

                << "Current Loss" << "," 
                << "Best Global Loss" << "," 
                << "Temperature" << "," 

                << "Accepted" << "," 
                << "Local improvement" << "," 
                << "Global improvement" << ","
                << "Stagnated" << "," 

                << "mutateStrokes time" << "," 
                << "render time" << "," 
                << "computeLoss time" << "," 
                << "minimizeColorError time" << "," 
 
                << "Iteration Time" << "\n";
            logf.close();
            std::cout << "OK: guardado log " << log_path << "\n";
        }
    }

    model.optimizeSimulatedAnnealing(1000000);

    output_filename = model.log_file_name + "_finalResult.png";
    if (!savePNG(model.getCurrentCanvas(), "output/images/" + output_filename)) {
        std::cerr << "Error guardando output/" << output_filename << "\n";
        return 1;
    }
    std::cout << "OK: guardado output/" << output_filename << "\n";

    return 0;
}
