#include <vector>
#include <string>
#include <cstdint>
#include <random>
#include "stroke.h"
#include <iostream>
#include <chrono>
#include <ctime>

// ================= Model ==================
// Class that manages the simulated annealing optimization model for the image aproximation problem.
// It has:
// - Variables: strokes vector, Current canvas, target image, for each stroke it saves a canvas size array with the stroke alphas drawn on it.
// - Methods: to render the current canvas, to compute the loss with respect to the target image, etc.
//              to start the initial greedy guess, to start the simulated annealing optimization, etc.
class Model {
private:
    std::vector<Stroke> strokes;
    Canvas currentCanvas;
    Canvas targetImage;
    ImageGray targetImageBorders;
    std::vector<ImageGray> brushes;

    std::vector<float> mutationWeights; // Weights for mutating x_rel, y_rel, size_rel, rotation_deg, color

    int n_strokes = 0;
    int n_sample_greedy = 1;
    float p = 0.1f;

    void minimizeColorError();
    void render_all();
    float computeLoss();

    std::string mutateStrokes();
public:
    Model(int n_strokes = 10, 
        const std::string& target_filename = "mona.png",
        int n_sample_greedy = 300,
        float p = 0.1f);
    ~Model();

    float max_loss = 0.0f;
    float current_loss = 0.0f;

    // LOGGING
    std::string target_name;
    std::string start_time;
    std::string last_time;

    std::string log_file_name;


    Canvas getCurrentCanvas() const { return currentCanvas; }

    void initialGuess();
    void optimizeSimulatedAnnealing(int n_iterations = 10000);
};