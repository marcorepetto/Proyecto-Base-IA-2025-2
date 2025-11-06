#include <vector>
#include <string>
#include <cstdint>
#include "stroke.h"

// ================= Model ==================
// Class that manages the optimization model for the image aproximation problem.
// It has:
// - Variables: strokes vector, Current canvas, target image, for each stroke it saves a canvas size array with the stroke alphas drawn on it.
// - Methods: to render the current canvas, to compute the loss with respect to the target image, etc.
// 
class Model {
private:
    std::vector<Stroke> strokes;
    Canvas currentCanvas;
    Canvas targetImage;
    ImageGray targetImageBorders;
    std::vector<ImageGray> brushes;
    float loss = 0.0f;

    int n_strokes = 0;
    int n_sample_greedy = 1;
    float pixel_threshold = 0.1f;
    float p = 0.1f;

    void minimizeColorError();
public:
    Model(int n_strokes = 10, 
        const std::string& target_filename = "mona.png",
        int n_sample_greedy = 300,
        float pixel_threshold = 0.1f,
        float p = 0.1f);
    ~Model();

    float max_loss = 0.0f;
    float current_loss = 0.0f;

    Canvas getCurrentCanvas() const { return currentCanvas; }
    
    void render();
    float computeLoss();
    void initialGuess();
};