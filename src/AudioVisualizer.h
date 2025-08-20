#pragma o*nce
#include "ofMain.h"
#include "AudioAnalyzer.h"
#include "Wave.h"
#include <memory>
#include <vector>
#include <map>

enum class VisualizationMode {
    WAVE_RMS,           // Raw waveform display
    SPECTRUM_BARS,      // FFT spectrum as bars
    CIRCULAR_SPECTRUM,  // FFT spectrum in circular form
    PEAK_PULSES,        // Visual pulses on peak detection
    CENTROID_WAVE,      // Wave colored by spectral centroid
    ONSET_PARTICLES,    // Particle system triggered by onsets
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

    // Clean separation: update processes all data, draw only renders
    void update(const AudioFeatures& features, const std::vector<float>& leftChannel = {},
        const std::vector<float>& rightChannel = {});
    void draw(); // No parameters needed - all data is processed internally

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

    // ALL processed visualization data (similar to smoothedRMS)
    // RMS data
    float smoothedRMS = 0.0f;
    float leftRMS = 0.0f;
    float rightRMS = 0.0f;

    // Spectrum data (processed and ready for rendering)
    std::vector<float> smoothedSpectrum;
    std::vector<float> spectrumBars;           // Processed spectrum for bar visualization
    std::vector<glm::vec2> circularSpectrum;  // Polar coordinates for circular spectrum

    // Peak and onset data
    float pulseIntensity = 0.0f;
    bool currentPeakDetected = false;
    bool currentOnsetDetected = false;

    // Centroid data
    float smoothedCentroid = 0.0f;
    float normalizedCentroid = 0.0f;  // Mapped to 0-1 range for color mapping

    // Particle system data
    std::vector<glm::vec2> particles;
    std::vector<float> particleLifetimes;
    std::vector<ofColor> particleColors;

    // Waveform data (processed and downsampled for rendering)
    std::vector<float> processedLeftChannel;
    std::vector<float> processedRightChannel;

    // Visual state
    float rotationAngle = 0.0f;
    ofColor currentColor;
    ofColor centroidColor;

    // History management
    static const int HISTORY_SIZE = 256;
    std::vector<float> spectrumHistory;

    // Update methods - all data processing happens here
    void updateRMSValues(const AudioFeatures& features);
    void updateSpectrumData(const AudioFeatures& features);
    void updatePeakData(const AudioFeatures& features);
    void updateOnsetData(const AudioFeatures& features);
    void updateCentroidData(const AudioFeatures& features);
    void updateParticles(bool addNew = false);
    void updateVisualState();
    void updateWaveformData(const std::vector<float>& leftChannel,
        const std::vector<float>& rightChannel);

    // Drawing methods - only rendering, no data processing
    void drawWaveformRMS();
    void drawSpectrumBars();
    void drawCircularSpectrum();
    void drawPeakPulses();
    void drawCentroidWave();
    void drawOnsetParticles();
    void drawCombinedView();

    // Helper functions
    ofColor getColorFromFrequency(float frequency);
    ofColor getColorFromAmplitude(float amplitude);
    void drawLabel(const std::string& text, glm::vec2 position);
    void addToHistory(float value, std::vector<float>& history);

    // Data processing helpers
    std::vector<float> downsampleForVisualization(const std::vector<float>& input, int targetSize);
    void smoothValue(float& target, float newValue, float smoothing);
    glm::vec2 polarToCartesian(float radius, float angle);
};