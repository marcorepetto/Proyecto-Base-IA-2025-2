#include "model.h"
#include "stroke.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <chrono>

#include <Eigen/Dense>


static std::mt19937 rng(std::random_device{}());
int randInt(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
}

// Model
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

    strokes[0].r = 1.0f; strokes[0].g = 0.0f; strokes[0].b = 0.0f; // rojo

    // Render inicial
    render();
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

    for (int i = 0; i < n_strokes; ++i) {
        Stroke& stroke = strokes[i];
        float best_score = 1.0f;
        Stroke best_stroke = stroke;
        int uncovered_pixels = 0;
        for (int j = 0; j < n_sample_greedy; ++j) {
            stroke.randomize();
            stroke.draw(currentCanvas);

            float score = 0.0f;
            int pixel_count = 0;
            int covered_pixels = 0;

            for (int y = 0; y < currentCanvas.height; ++y) {
                for (int x = 0; x < currentCanvas.width; ++x) {
                    float a = stroke.strokeAlphas[y * currentCanvas.width + x];
                    if (a <= 0.0f) {
                        if (targetImageBorders.gray[y * currentCanvas.width + x] < pixel_threshold) {
                            uncovered_pixels++;
                        }
                        continue;
                    }
                    pixel_count++;
                    int idx = y * currentCanvas.width + x;
                    covered_pixels += targetImageBorders.gray[idx] > 0.0f ? 1 : 0;
                }
            }

            float uncovered_ratio = float(uncovered_pixels) / float(n_pixels - pixel_count);
            float already_covered_ratio = float(covered_pixels) / float(pixel_count);
            
            score = (1.0f - p) * uncovered_ratio + p * already_covered_ratio;

            if (score < best_score) {
                best_score = score;
                best_stroke = stroke;

                std::cout << "["<< i+1 << "] ["<< j+1 << "] New best stroke [Score: " << best_score << "]: "
                        << "x_rel=" << best_stroke.x_rel << ", y_rel=" << best_stroke.y_rel
                        << ", size_rel=" << best_stroke.size_rel << ", rotation_deg=" << best_stroke.rotation_deg
                        << ", type=" << best_stroke.type << "\n";
                
                continue;
            }
            std::cout << "["<< i+1<< "] ["<< j+1 << "]                 [Score: " << best_score << "]\n";
        }
        stroke = best_stroke;
        for (int y = 0; y < currentCanvas.height; ++y) {
            for (int x = 0; x < currentCanvas.width; ++x) {
                int idx = y * currentCanvas.width + x;
                targetImageBorders.gray[idx] = stroke.strokeAlphas[idx] > 0.0f ? pixel_threshold : targetImageBorders.gray[idx];
            }
        }
    }

    minimizeColorError();
}

