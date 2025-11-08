#include "model.h"
#include "stroke.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include <limits>

#include <Eigen/Dense>


static std::mt19937 rng(std::random_device{}());
int randInt(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
}

/**
 * @brief Constructor de la clase Model.
 * @param n_strokes Número de strokes a utilizar.
 * @param target_filename Nombre del archivo de la imagen objetivo.
 * @param n_sample_greedy Número de muestras aleatorias por stroke para la inicialización
 * @param pixel_threshold Umbral de borde para la inicialización.
 * @param p Peso del borde respecto al área cubierta en la inicialización.
 */
Model::Model(int n_strokes, 
            const std::string& target_filename,
            int n_sample_greedy,
            float pixel_threshold,
            float p) : 
        currentCanvas(1, 1), 
        targetImage(1, 1), 
        n_strokes(n_strokes),
        n_sample_greedy(n_sample_greedy),
        pixel_threshold(pixel_threshold),
        p(p) {

    // Cargar imagen objetivo (target)
    std::cout << "Cargando imagen objetivo " << target_filename << "...\n";
    if (!loadImageRGB_asCanvas("instancias/" + target_filename, targetImage)) {
        std::cerr << "Model: no pude cargar " << target_filename << "\n";
    }

    const int canvas_width = targetImage.width;
    const int canvas_height = targetImage.height;
    
    // Target image borders
    loadImageGray("borders/" + target_filename, targetImageBorders);

    // multiplicar por 1/255.0f
    for (size_t i = 0; i < targetImageBorders.gray.size(); ++i) {
        targetImageBorders.gray[i] *= 0.1f;
    }

    currentCanvas = Canvas(canvas_width, canvas_height);

    // Calculate max_loss
    max_loss = 0.0f;
    const int N = canvas_width * canvas_height * 3;
    for (int i = 0; i < N; ++i) {
        float pixel = targetImage.rgb[i];
        float max_pixel_error = std::max(pixel * pixel, (1.0f - pixel) * (1.0f - pixel));
        max_loss += max_pixel_error;
    }
    max_loss = max_loss / (canvas_width * canvas_height);

    // Inicializar strokes
    strokes.reserve(n_strokes);
    for (int i = 0; i < n_strokes; ++i) {
        float xr = 0.0f;
        float yr = 0.0f;
        float sr = 0.3f; // tamaño relativo
        float rot = 0;
        int type = 0;
        float rr = 1.0f;
        float gg = 1.0f;
        float bb = 1.0f;

        Stroke S(xr, yr, sr, rot, type, rr, gg, bb, canvas_width, canvas_height);
        S.draw(currentCanvas);
        strokes.push_back(S);
    }
}

Model::~Model() {}

float Model::computeLoss() {
    // Calcular la pérdida (loss) entre el canvas actual y la imagen objetivo
    float loss = 0.0f;
    const int N = currentCanvas.width * currentCanvas.height * 3;
    for (int i = 0; i < N; ++i) {
        float diff = float(currentCanvas.rgb[i]) - float(targetImage.rgb[i]);
        loss += diff * diff;
    }
    loss = loss / (currentCanvas.width * currentCanvas.height);
    return loss;
}

void Model::render() {
    currentCanvas.clear(0, 0, 0);
    Stroke stroke;
    
    for (int y = 0; y < currentCanvas.height; ++y) {
        for (int x = 0; x < currentCanvas.width; ++x) {
            float a_factor = 1.0f;
            for (int i = (int)strokes.size() - 1; i >= 0; --i) {
                stroke = strokes[i];
                float a = stroke.strokeAlphas[y * currentCanvas.width + x];
                if (a <= 0.0f) continue;

                // Mezcla alfa simple
                int idx = (y * currentCanvas.width + x) * 3;

                float r_weight = (stroke.r * a * a_factor);
                float g_weight = (stroke.g * a * a_factor);
                float b_weight = (stroke.b * a * a_factor);

                currentCanvas.rgb[idx + 0] += r_weight;
                currentCanvas.rgb[idx + 1] += g_weight;
                currentCanvas.rgb[idx + 2] += b_weight;

                a_factor *= (1.0f - a);
                if (a_factor <= 0.001f) break; // casi opaco, salir
            }
        }
    }
}

void Model::minimizeColorError() {
    // Minimizar error **POR COLOR** para cada stroke en conjunto mediante minimos cuadrados
    Eigen::MatrixXf A(currentCanvas.width * currentCanvas.height, n_strokes);
    Eigen::VectorXf b_r(currentCanvas.width * currentCanvas.height);
    Eigen::VectorXf b_g(currentCanvas.width * currentCanvas.height);
    Eigen::VectorXf b_b(currentCanvas.width * currentCanvas.height);
    A.setZero();
    b_r.setZero();
    b_g.setZero();
    b_b.setZero();

    std::cout << "Construyendo sistema de ecuaciones para mínimos cuadrados...\n";
    for (int y = 0; y < currentCanvas.height; ++y) {
        std::cout << "  Procesando fila " << y+1 << " de " << currentCanvas.height << "...\n";
        for (int x = 0; x < currentCanvas.width; ++x) {
            // Fill A
            float alpha_accum = 1.0f;
            for (int i = (int)strokes.size() - 1; i >= 0; --i) {
                const Stroke& stroke = strokes[i];
                float a = stroke.strokeAlphas[y * currentCanvas.width + x];
                A(y * currentCanvas.width + x, i) = a * alpha_accum;
                alpha_accum *= (1.0f - a);
            }

            // Fill b
            int idx = (y * currentCanvas.width + x) * 3;
            b_r(y * currentCanvas.width + x) = targetImage.rgb[idx + 0];
            b_g(y * currentCanvas.width + x) = targetImage.rgb[idx + 1];
            b_b(y * currentCanvas.width + x) = targetImage.rgb[idx + 2];
        }
    }
    std::cout << "Sistema construido: A es " << A.rows() << "x" << A.cols() << "\n\nRealizando descomposición QR de A (" << A.rows() << "x" << A.cols() << ")... ";

    // QR decomposition to solve for r, g, b. Compute QR = A once and reuse.
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::ColPivHouseholderQR<Eigen::MatrixXf> qr(A);

    auto end = std::chrono::high_resolution_clock::now();
    // Save seconds duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "hecho en " << duration / 1000000.0 << " s.\n";

    std::cout << "Resolviendo para canales R, G, B...\n";
    start = std::chrono::high_resolution_clock::now();
    std::cout << "  Canal R... ";
    Eigen::VectorXf x_r = qr.solve(b_r);
    std::cout << "Listo!\n  Canal G... ";
    Eigen::VectorXf x_g = qr.solve(b_g);
    std::cout << "Listo!\n  Canal B... ";
    Eigen::VectorXf x_b = qr.solve(b_b);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Listo!\nHecho en " << duration / 1000000.0 << " s.\n";

    /*
    std::cout << "Matrix A:\n" << A << "\n";
    std::cout << "Vector b_r:\n" << b_r << "\n";
    std::cout << "Vector b_g:\n" << b_g << "\n";
    std::cout << "Vector b_b:\n" << b_b << "\n";

    std::cout << "Initial guess colors computed via least squares.\n";
    std::cout << "Red channel approximation using x_r (then real values):\n" 
        << x_r.transpose() << "\n" 
        << (A*x_r).transpose() << "\n" 
        << b_r.transpose() << "\n";
    std::cout << "Green channel approximation using x_g (then real values):\n" 
        << x_g.transpose() << "\n" 
        << b_g.transpose() << "\n";
    std::cout << "Blue channel approximation using x_b (then real values):\n" 
        << x_b.transpose() << "\n" 
        << (A*x_b).transpose() << "\n" 
        << b_b.transpose() << "\n";

    std::cout << "Colors:\n" 
        << x_r.transpose() << "\n" 
        << x_g.transpose() << "\n" 
        << x_b.transpose() << "\n";
    */

    // Update stroke colors
    for (int i = 0; i < n_strokes; ++i) {
        strokes[i].r = std::clamp(x_r(i), 0.0f, 1.0f);
        strokes[i].g = std::clamp(x_g(i), 0.0f, 1.0f);
        strokes[i].b = std::clamp(x_b(i), 0.0f, 1.0f);
    }
}

void Model::initialGuess() {
    int n_pixels = currentCanvas.width * currentCanvas.height;
    const float alpha_structural = 0.98f; // Weight for gradient
    const float beta_texture = 0.05f;     // Weight for texture

    // 1. Create structural energy map E(x,y)
    // 1.1 Convert targetImage to grayscale
    ImageGray targetGray;
    targetGray.width = targetImage.width;
    targetGray.height = targetImage.height;
    targetGray.gray.resize(n_pixels);
    for (int i = 0; i < n_pixels; ++i) {
        float r = targetImage.rgb[i*3 + 0];
        float g = targetImage.rgb[i*3 + 1];
        float b = targetImage.rgb[i*3 + 2];
        targetGray.gray[i] = 0.299f * r + 0.587f * g + 0.114f * b;
    }

    // 1.2 Compute texture map T(x,y) = local variance
    ImageGray textureMap;
    textureMap.width = targetImage.width;
    textureMap.height = targetImage.height;
    textureMap.gray.resize(n_pixels);
    const int w_size = 5;
    const int w_half = w_size / 2;
    for (int y = 0; y < targetGray.height; ++y) {
        for (int x = 0; x < targetGray.width; ++x) {
            float sum = 0;
            float sum_sq = 0;
            int count = 0;
            for (int wy = -w_half; wy <= w_half; ++wy) {
                for (int wx = -w_half; wx <= w_half; ++wx) {
                    int cx = x + wx;
                    int cy = y + wy;
                    if (cx >= 0 && cx < targetGray.width && cy >= 0 && cy < targetGray.height) {
                        float val = targetGray.gray[cy * targetGray.width + cx];
                        sum += val;
                        sum_sq += val * val;
                        count++;
                    }
                }
            }
            if (count > 0) {
                float mean = sum / count;
                float variance = (sum_sq / count) - (mean * mean);
                textureMap.gray[y * targetGray.width + x] = std::sqrt(std::max(0.0f, variance));
            } else {
                textureMap.gray[y * targetGray.width + x] = 0.0f;
            }
        }
    }
    
    // Normalize texture map
    float max_texture = 0.0f;
    for(float val : textureMap.gray) max_texture = std::max(max_texture, val);
    if (max_texture > 0) {
        for(float& val : textureMap.gray) val /= max_texture;
    }

    // Normalize gradient map
    ImageGray gradientMap = targetImageBorders; // Copy
    float max_gradient = 0.0f;
    for(float val : gradientMap.gray) max_gradient = std::max(max_gradient, val);
    if (max_gradient > 0) {
        for(float& val : gradientMap.gray) val /= max_gradient;
    }

    // 1.3 Combine G and T to get E
    ImageGray energyMap;
    energyMap.width = targetImage.width;
    energyMap.height = targetImage.height;
    energyMap.gray.resize(n_pixels);

    for (int i = 0; i < n_pixels; ++i) {
        energyMap.gray[i] = alpha_structural * gradientMap.gray[i] + beta_texture * textureMap.gray[i];
    }
    // Normalize energyMap
    float max_energy = 0.0f;
    for(float val : energyMap.gray) max_energy = std::max(max_energy, val);
    if (max_energy > 0) {
        for(float& val : energyMap.gray) val /= max_energy;
    }

    // Coverage map C(p)
    std::vector<float> coverageMap(n_pixels, 0.0f);

    // Main loop
    for (int i = 0; i < n_strokes; ++i) {
        Stroke& stroke = strokes[i];
        float best_score = -std::numeric_limits<float>::max();
        Stroke best_stroke = stroke;

        for (int j = 0; j < n_sample_greedy; ++j) {
            stroke.randomize();
            
            stroke.draw(currentCanvas); // This is needed to compute strokeAlphas

            float gain = 0.0f;
            float penalty = 0.0f;
            float sum_a = 0.0f;

            for (int p_idx = 0; p_idx < n_pixels; ++p_idx) {
                float a = stroke.strokeAlphas[p_idx];
                sum_a += a;
                if (a > 0.0f) {
                    gain += a * energyMap.gray[p_idx];
                    penalty += a * coverageMap[p_idx];
                }
            }

            float score = gain/sum_a - p * penalty; // Use p as lambda

            if (score > best_score) {
                best_score = score;
                best_stroke = stroke;

                if (j % 10 == 0) { // Log less frequently
                    std::cout << "["<< i+1 << "] ["<< j+1 << "] New best stroke [Score: " << best_score << "]: "
                            << "gain = " << gain << ", "
                            << "penalty = " << penalty << "\n";
                }
            }
        }
        stroke = best_stroke;
        float scale_factor = 0.1f;
        stroke.size_rel = stroke.size_rel + scale_factor > 0.5f ? 0.5f : stroke.size_rel + scale_factor; // Increase size a bit for better coverage
        stroke.draw(currentCanvas); // Ensure best_stroke alphas are calculated

        // Update energy map and coverage map
        for (int p_idx = 0; p_idx < n_pixels; ++p_idx) {
            float a = stroke.strokeAlphas[p_idx];
            if (a > 0.0f) {
                energyMap.gray[p_idx] *= (1.0f - a);
                coverageMap[p_idx] += a;
            }
        }
    }

    minimizeColorError();
    render();
}

void Model::optimizeSimulatedAnnealing(int n_iterations) {
    // Simulated annealing variables
    float current_loss = computeLoss();
    std::vector<Stroke> prev_strokes;
    Canvas prev_canvas = currentCanvas;

    // Initial temperature defined by max_loss
    float T_initial = max_loss;
    float T_final = 0.1f;
    float T = T_initial;
    float alpha = std::pow(T_final / T_initial, 1.0f / n_iterations);

    // Restart variables
    float best_global_loss = current_loss;
    std::vector<Stroke> best_global_strokes = strokes;
    Canvas best_global_canvas = currentCanvas;
    int stagnation_counter = 0;
    const int stagnation_limit = n_iterations / 10;

    // Random engine
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<float> U(0.0f, 1.0f);

    for (int iter = 0; iter < n_iterations; ++iter) {
        // Save current state
        prev_strokes = strokes;
        prev_canvas = currentCanvas;

        // Mutate strokes
        mutateStrokes();

        // Render new canvas
        render();

        // Compute new loss
        float new_loss = computeLoss();

        // Acceptance probability
        float delta_loss = new_loss - current_loss;
        if (delta_loss < 0 || std::exp(-delta_loss / T) > U(gen)) {
            // Accept new state
            current_loss = new_loss;

            // Update best global solution
            if (current_loss < best_global_loss) {
                best_global_loss = current_loss;
                best_global_strokes = strokes;
                best_global_canvas = currentCanvas;
                stagnation_counter = 0;
            } else {
                stagnation_counter++;
            }
        } else {
            // Revert to previous state
            strokes = prev_strokes;
            currentCanvas = prev_canvas;
        }

        // Update temperature
        T *= alpha;

        // Check for stagnation
        if (stagnation_counter >= stagnation_limit) {
            std::cout << "Stagnation detected at iteration " << iter << ". Restarting from best solution.\n";
            strokes = best_global_strokes;
            currentCanvas = best_global_canvas;
            current_loss = best_global_loss;
            stagnation_counter = 0;
        }

        // Logging
        if (iter % 100 == 0) {
            std::cout << "Iteration " << iter << ": Current Loss = " << current_loss 
                      << ", Best Global Loss = " << best_global_loss << ", Temperature = " << T << "\n";
        }
    }
}

void Model::mutateStrokes() {
    
}