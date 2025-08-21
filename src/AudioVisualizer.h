#pragma once
#include "ofMain.h"
#include "AudioAnalyzer.h"
#include "Wave.h"
#include <memory>
#include <vector>
#include <map>

struct VisualizationParams {
    glm::vec2 position = { 0, 0 };
    glm::vec2 size = { 400, 300 };
    ofColor primaryColor = ofColor::white;
    ofColor secondaryColor = ofColor::gray;
    float scale = 1.0f;
    float smoothing = 0.0001f;  // Note: This might be unused now
    bool showLabels = true;
    float sensitivity = 1.0f;
};

class AudioVisualizer {
public:
    AudioVisualizer();
    ~AudioVisualizer();

    void setup(int sampleRate, int bufferSize);
    void setMode(VisualizationMode mode);
    void setParams(const VisualizationParams& params) { visualParams = params; }

    void update(const AudioFeatures& features);

    void draw(); // Clean draw method - no parameters needed

    // Mode management
    VisualizationMode getCurrentMode() const { return currentMode; }
    std::vector<std::string> getModeNames() const;
    void nextMode();
    void previousMode();

    // Parameter control
    VisualizationParams& getParams() { return visualParams; }
    const VisualizationParams& getParams() const { return visualParams; }

    // Utility functions
    void reset();

private:
    VisualizationMode currentMode;
    VisualizationParams visualParams;

    int sampleRate;
    int bufferSize;

    // Current audio features (all processed data comes from here)
    AudioFeatures currentFeatures;

    // Visualization-specific components
    std::vector<Wave> waves;  // For centroid wave visualization

    // Particle system (only remaining local state that needs updating)
    std::vector<glm::vec2> particles;
    std::vector<float> particleLifetimes;
    std::vector<ofColor> particleColors;

    // Drawing methods - only rendering, no data processing
    void drawWaveformRMS();
    void drawSpectrumBars();
    void drawCircularSpectrum();
    void drawPeakPulses();
    void drawCentroidWave();
    void drawOnsetParticles();
    void drawCombinedView();

    // Particle management (only remaining update logic)
    void updateParticles(bool addNew);

    // Helper functions for rendering
    ofColor getColorFromFrequency(float frequency);
    ofColor getColorFromAmplitude(float amplitude);
};