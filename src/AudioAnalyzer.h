// AudioAnalyzer.h
#pragma once
#include "ofMain.h"
#include <memory>
#include <vector>
#include <complex>


enum class AudioAnalysisMode {
    RMS_AMPLITUDE,
    FFT_SPECTRUM,
    PEAK_DETECTION,
    SPECTRAL_CENTROID,
    ONSET_DETECTION
};

struct AudioFeatures {
    // RMS data
    float leftRMS = 0.0f;
    float rightRMS = 0.0f;

    // FFT data
    std::vector<float> fftMagnitudes;
    std::vector<float> fftBins;

    // Peak detection
    bool peakDetected = false;
    float peakMagnitude = 0.0f;

    // Spectral centroid
    float spectralCentroid = 0.0f;

    // Onset detection
    bool onsetDetected = false;
    float onsetStrength = 0.0f;

    // General
    float overallEnergy = 0.0f;};

class AudioAnalyzer {
public:
    AudioAnalyzer();
    ~AudioAnalyzer();

    void setup(int sampleRate, int bufferSize);
    void setMode(AudioAnalysisMode mode);
    AudioAnalysisMode getCurrentMode() const { return currentMode; }

    // Main analysis function
    AudioFeatures analyze(const std::vector<float>& leftChannel,
        const std::vector<float>& rightChannel);

    // Mode-specific analysis functions
    AudioFeatures analyzeRMS(const std::vector<float>& left, const std::vector<float>& right);
    AudioFeatures analyzeFFT(const std::vector<float>& left, const std::vector<float>& right);
    AudioFeatures analyzePeaks(const std::vector<float>& left, const std::vector<float>& right);
    AudioFeatures analyzeSpectralCentroid(const std::vector<float>& left, const std::vector<float>& right);
    AudioFeatures analyzeOnsets(const std::vector<float>& left, const std::vector<float>& right);

    // Utility functions
    std::vector<std::string> getModeNames() const;
    void nextMode();
    void previousMode();

    void setAutoScaling(bool enable) { enableAutoScaling = enable; }
    float getMaxObservedRMS() const { return maxObservedRMS; }
    void resetLevelHistory() { maxObservedRMS = 0.001f; }


private:
    AudioAnalysisMode currentMode;
    int sampleRate;
    int bufferSize;
    float maxObservedRMS = 0.001f; // Start with small value to avoid division by zero
    float adaptationRate = 0.001f; // How fast to adapt to new levels
    bool enableAutoScaling = true;

    // FFT-related members
    std::vector<std::complex<float>> fftInput;
    std::vector<std::complex<float>> fftOutput;
    void performFFT(const std::vector<float>& input, std::vector<float>& magnitudes);

    // Peak detection members
    std::vector<float> previousSamples;
    float peakThreshold;

    // Onset detection members
    std::vector<float> previousSpectrum;
    float onsetThreshold;

    // Smoothing for various analyses
    float smoothingFactor;
    AudioFeatures previousFeatures;

    // Helper functions
    float calculateRMS(const std::vector<float>& samples);
    float calculateSpectralCentroid(const std::vector<float>& magnitudes, const std::vector<float>& frequencies);
    bool detectPeak(const std::vector<float>& samples, float threshold);
    bool detectOnset(const std::vector<float>& currentSpectrum, const std::vector<float>& previousSpectrum);
};