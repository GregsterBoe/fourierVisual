// AudioAnalyzer.h
#pragma once
#include "ofMain.h"
#include <memory>
#include <vector>
#include <complex>
#include <cmath>


// Simplified structures - remove rhythm/accent complexity
struct MelodyRange {
    float start;
    float end;
    float centroid;
    float energy;
    float confidence;

    MelodyRange() : start(0), end(0), centroid(0), energy(0), confidence(0) {}
    MelodyRange(float s, float e, float c, float en) : start(s), end(e), centroid(c), energy(en), confidence(0) {}
};

class MelodyTracker {
private:
    // Configuration for note detection
    float adaptationRate = 0.25f;        // Fast adaptation for note changes
    float energyThreshold = 0.02f;       // Lower threshold for subtle notes
    float confidenceDecay = 0.92f;       // Faster decay to track note changes
    float minConfidence = 0.2f;          // Lower confidence threshold

    // Melody frequency bounds (where most musical notes occur)
    float melodyMin = 80.0f;             // Around low E2 (guitar)
    float melodyMax = 4000.0f;           // Above high soprano range

    // State
    MelodyRange currentRange;
    std::vector<float> energyHistory;
    std::vector<float> previousBands;    // For smoothing
    int historyFrames = 10;              // Shorter history for note tracking

    // Helper methods
    float calculateMelodyCentroid(const std::vector<float>& magnitudes, float binWidth);
    float calculateMelodyEnergy(const std::vector<float>& magnitudes, float binWidth);
    void updateMelodyRange(const std::vector<float>& magnitudes, float binWidth);
    std::pair<float, float> findMelodyPeaks(const std::vector<float>& magnitudes, float binWidth);

public:
    MelodyTracker() {
        energyHistory.reserve(historyFrames);
        reset();
    }

    void reset();
    void updateRange(const std::vector<float>& magnitudes, float sampleRate);
    std::vector<float> getMelodyBands(const std::vector<float>& magnitudes, float sampleRate, int numBands);
    std::pair<float, float> getCurrentMelodyRange() const;
    float getDominantFrequency() const { return currentRange.centroid; }
    float getMelodyConfidence() const { return currentRange.confidence; }

    void setAdaptationRate(float rate) { adaptationRate = std::max(0.01f, std::min(1.0f, rate)); }
};

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

    std::vector<float> logFrequencyBands;   // Logarithmically spaced bands
    float bassEnergy = 0.0f;               // Overall bass energy
    float midEnergy = 0.0f;                // Mid-range energy
    float trebleEnergy = 0.0f;             // High-frequency energy


    // General
    float overallEnergy = 0.0f;

    // Simplified melody tracking
    float dominantFrequency = 0.0f;     // The detected note frequency
    float melodyConfidence = 0.0f;      // How confident we are about the note
    std::pair<float, float> melodyRange = { 0.0f, 0.0f }; // Current tracking range
};

class AudioAnalyzer {
public:
    AudioAnalyzer();
    ~AudioAnalyzer();
    bool enableMelodyTracking = false;


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


    void resetMelodyTracking() {
        melodyTracker.reset();
    }

    std::pair<float, float> getCurrentMelodyRange() const {
        return melodyTracker.getCurrentMelodyRange();
    }

    float getDominantMelodyFrequency() const {
        return melodyTracker.getDominantFrequency();
    }

    float getMelodyConfidence() const {
        return melodyTracker.getMelodyConfidence();
    }

    void setMelodyAdaptationRate(float rate) {
        melodyTracker.setAdaptationRate(rate);
    }

private:
    // Basic settings
    int sampleRate;
    int bufferSize;

    // Auto-scaling parameters
    float maxObservedRMS = 0.001f;
    float adaptationRate = 0.001f;
    bool enableAutoScaling = true;
    float runningMax = 0.0f;
    float maxDecayRate = 0.99f;

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
    void performFFTFast(const std::vector<float>& input, std::vector<float>& magnitudes);
    void performRealFFT(const std::vector<float>& input, std::vector<float>& magnitudes);
    void performFFTOptimal(const std::vector<float>& input, std::vector<float>& magnitudes);
    void cooleyTukeyFFT(std::vector<std::complex<float>>& data);

    // Frequency response enhancement methods
    void enhanceFrequencySelectivity(std::vector<float>& magnitudes);

    void applySelectiveSmoothing(std::vector<float>& magnitudes);


    // Enhanced frequency analysis
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

    MelodyTracker melodyTracker;
};