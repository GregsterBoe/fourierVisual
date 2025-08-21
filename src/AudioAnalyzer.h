// AudioAnalyzer.h
#pragma once
#include "ofMain.h"
#include <memory>
#include <vector>
#include <complex>
#include <cmath>

enum class VisualizationMode {
    WAVE_RMS,           // Raw waveform display
    SPECTRUM_BARS,      // FFT spectrum as bars
    CIRCULAR_SPECTRUM,  // FFT spectrum in circular form
    PEAK_PULSES,        // Visual pulses on peak detection
    CENTROID_WAVE,      // Wave colored by spectral centroid
    ONSET_PARTICLES,    // Particle system triggered by onsets
    COMBINED_VIEW       // Multiple visualizations combined
};

struct AudioFeatures {
    // RMS data
    float leftRMS = 0.0f;
    float rightRMS = 0.0f;
    float smoothedRMS = 0.0f;  // Pre-smoothed for visualization

    // FFT data - raw and processed
    std::vector<float> fftMagnitudes;
    std::vector<float> fftBins;
    std::vector<float> smoothedSpectrum;     // Pre-smoothed spectrum
    std::vector<float> spectrumBarHeights;   // Ready-to-render bar heights
    std::vector<glm::vec2> circularSpectrum; // Pre-calculated polar coordinates

    // Peak detection
    bool peakDetected = false;
    float peakMagnitude = 0.0f;
    float pulseIntensity = 0.0f;  // Pre-calculated pulse intensity with decay

    // Spectral centroid
    float spectralCentroid = 0.0f;
    float smoothedCentroid = 0.0f;     // Pre-smoothed
    float normalizedCentroid = 0.0f;   // Normalized 0-1 for color mapping
    ofColor centroidColor;             // Pre-calculated color

    // Onset detection
    bool onsetDetected = false;
    float onsetStrength = 0.0f;

    // Waveform data (downsampled for visualization)
    std::vector<float> leftChannelViz;   // Downsampled left channel
    std::vector<float> rightChannelViz;  // Downsampled right channel

    // Visual state
    float rotationAngle = 0.0f;     // For rotating visualizations
    ofColor currentColor;           // Pre-calculated color based on audio

    std::vector<float> bassResponse;        // Detailed low-frequency bands
    std::vector<float> logFrequencyBands;   // Logarithmically spaced bands
    float bassEnergy = 0.0f;               // Overall bass energy
    float midEnergy = 0.0f;                // Mid-range energy
    float trebleEnergy = 0.0f;             // High-frequency energy

    // General
    float overallEnergy = 0.0f;
};

class AudioAnalyzer {
public:
    AudioAnalyzer();
    ~AudioAnalyzer();

    void setup(int sampleRate, int bufferSize);

    // Main analysis function - now provides comprehensive features
    AudioFeatures analyze(const std::vector<float>& leftChannel,
        const std::vector<float>& rightChannel);

    // Auto-scaling controls
    void setAutoScaling(bool enable) { enableAutoScaling = enable; }
    float getMaxObservedRMS() const { return maxObservedRMS; }
    float getMaxFrequency() const {
        return maxFrequency;
    }
    float getAverageFrequency() const {
        return averageFrequency;
    }

    void resetLevelHistory() { maxObservedRMS = 0.001f; }

    // Visualization synchronization
    void setVisualizationSize(glm::vec2 size) { visualizationSize = size; }

private:
    // Basic settings
    int sampleRate;
    int bufferSize;

    // Auto-scaling parameters
    float maxObservedRMS = 0.001f;
    float adaptationRate = 0.001f;
    bool enableAutoScaling = true;

    // frequency limits
    float lowFrequencyLimit = 250.0f;
    float midFrequencyLimit = 4000.0f;
    float highFrequencyLimit = 16000.0f;;

    float lowFrequencyScaling = 5.0f;

    float frequencySum = 0.0f;
    float maxFrequency = 0.01f;
    float averageFrequency = 0.01f;

    // Analysis thresholds
    float peakThreshold;
    float onsetThreshold;
    float smoothingFactor;

    // FFT buffers
    std::vector<std::complex<float>> fftInput;
    std::vector<std::complex<float>> fftOutput;

    // State tracking for analysis
    std::vector<float> previousSamples;
    std::vector<float> previousSpectrum;
    std::vector<float> previousMagnitudes;

    // Visualization processing state
    float smoothedRMS = 0.0f;
    float smoothedCentroid = 0.0f;
    float pulseIntensity = 0.0f;
    float rotationAngle = 0.0f;
    std::vector<float> smoothedSpectrum;

    // Visualization parameters
    float visualSmoothingFactor = 0.05f;
    glm::vec2 visualizationSize = { 400, 300 };

    // Core analysis functions
    float calculateRMS(const std::vector<float>& samples);
    
    // fft methods 
    void performFFT(const std::vector<float>& input, std::vector<float>& magnitudes);
    void performFFTFast(const std::vector<float>& input, std::vector<float>& magnitudes);
    void performRealFFT(const std::vector<float>& input, std::vector<float>& magnitudes);
    void performFFTOptimal(const std::vector<float>& input, std::vector<float>& magnitudes);
    void cooleyTukeyFFT(std::vector<std::complex<float>>& data);

    // Frequency response enhancement methods
    float applyFrequencyCompensation(float magnitude, float frequency);
    float amplitudeToVisualizationScale(float magnitude, float frequency);
    void applyFrequencyAwareSmoothing(std::vector<float>& magnitudes);

    // Enhanced frequency analysis
    std::vector<float> getLowFrequencyBands(const std::vector<float>& magnitudes, int numBands = 4);
    std::vector<float> getLogFrequencyBands(const std::vector<float>& magnitudes, int numBands);
    void calculateFrequencyEnergies(AudioFeatures& features);

    float calculateSpectralCentroid(const std::vector<float>& magnitudes,
        const std::vector<float>& frequencies);
    bool detectPeak(const std::vector<float>& samples, float threshold);
    bool detectOnset(const std::vector<float>& currentSpectrum,
        const std::vector<float>& previousSpectrum);

    // Visualization processing helpers
    std::vector<float> downsampleForVisualization(const std::vector<float>& input, int targetSize);
    ofColor getColorFromFrequency(float freq);
    ofColor getColorFromAmplitude(float amplitude);
};