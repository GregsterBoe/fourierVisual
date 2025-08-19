#pragma once
#include "ofMain.h"
#include "AudioAnalyzer.h"
#include "Wave.h"
#include <memory>
#include <vector>
#include <map>

enum class VisualizationMode {
    WAVE_RMS,           // Classic wave with RMS amplitude
    SPECTRUM_BARS,      // FFT spectrum as bars
    CIRCULAR_SPECTRUM,  // FFT spectrum in circular form
    PEAK_PULSES,        // Visual pulses on peak detection
    CENTROID_WAVE,      // Wave colored by spectral centroid
    ONSET_PARTICLES,    // Particle system triggered by onsets
    WAVEFORM_RAW,       // Raw waveform display
    COMBINED_VIEW       // Multiple visualizations combined
};

struct VisualizationParams {
    glm::vec2 position = { 0, 0 };
    glm::vec2 size = { 400, 300 };
    ofColor primaryColor = ofColor::white;
    ofColor secondaryColor = ofColor::gray;
    float scale = 1.0f;
    float smoothing = 0.0001f;
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
    void update(const AudioFeatures& features, const std::vector<float>& leftChannel = {},
        const std::vector<float>& rightChannel = {});
    void draw(AudioFeatures& features);

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
    void setAutoScale(bool enable) { autoScale = enable; }

private:
    VisualizationMode currentMode;
    VisualizationParams visualParams;

    int sampleRate;
    int bufferSize;
    bool autoScale = true;
    int defaultScale = 20;

    // Visualization-specific data
    std::vector<Wave> waves;
    std::vector<float> spectrumHistory;
    std::vector<glm::vec2> particles;
    std::vector<float> particleLifetimes;

    // Smoothed values for visual stability
    float smoothedRMS = 0.0f;
    float leftRMS = 0.0f;
    float rightRMS = 0.0f;
    float smoothedCentroid = 0.0f;
    std::vector<float> smoothedSpectrum;

    // Visual state
    float pulseIntensity = 0.0f;
    float rotationAngle = 0.0f;
    ofColor currentColor;

    // Mode-specific drawing functions
    void drawWaveRMS();
    void drawSpectrumBars(const AudioFeatures& features);
    void drawCircularSpectrum(const AudioFeatures& features);
    void drawPeakPulses(const AudioFeatures& features);
    void drawCentroidWave(const AudioFeatures& features, const std::vector<float>& leftChannel);
    void drawOnsetParticles(const AudioFeatures& features);
    void drawRawWaveform(const std::vector<float>& leftChannel, const std::vector<float>& rightChannel);
    void drawCombinedView(const AudioFeatures& features, const std::vector<float>& leftChannel, const std::vector<float>& rightChannel);

    // Helper functions
    void updateParticles(bool addNew = false);
    void updateSmoothing(const AudioFeatures& features);
    ofColor getColorFromFrequency(float frequency);
    ofColor getColorFromAmplitude(float amplitude);
    void drawLabel(const std::string& text, glm::vec2 position);

    // History management for certain visualizations
    static const int HISTORY_SIZE = 256;
    void addToHistory(float value, std::vector<float>& history);
};