#include "AudioAnalyzer.h"
#include <algorithm>
#include <numeric>

AudioAnalyzer::AudioAnalyzer()
    : currentMode(AudioAnalysisMode::RMS_AMPLITUDE)
    , sampleRate(44100)
    , bufferSize(512)
    , peakThreshold(0.1f)
    , onsetThreshold(0.05f)
    , smoothingFactor(0.1f) {
}

AudioAnalyzer::~AudioAnalyzer() {
}

void AudioAnalyzer::setup(int sampleRate, int bufferSize) {
    this->sampleRate = sampleRate;
    this->bufferSize = bufferSize;

    // Initialize FFT buffers
    fftInput.resize(bufferSize);
    fftOutput.resize(bufferSize);

    // Initialize other buffers
    previousSamples.resize(bufferSize, 0.0f);
    previousSpectrum.resize(bufferSize / 2, 0.0f);
}

void AudioAnalyzer::setMode(AudioAnalysisMode mode) {
    currentMode = mode;
    ofLogNotice("AudioAnalyzer") << "Switched to mode: " << static_cast<int>(mode);
}

AudioFeatures AudioAnalyzer::analyze(const std::vector<float>& leftChannel,
    const std::vector<float>& rightChannel) {

    // Dispatch to appropriate analysis function based on current mode
    switch (currentMode) {
    case AudioAnalysisMode::RMS_AMPLITUDE:
        return analyzeRMS(leftChannel, rightChannel);

    case AudioAnalysisMode::FFT_SPECTRUM:
        return analyzeFFT(leftChannel, rightChannel);

    case AudioAnalysisMode::PEAK_DETECTION:
        return analyzePeaks(leftChannel, rightChannel);

    case AudioAnalysisMode::SPECTRAL_CENTROID:
        return analyzeSpectralCentroid(leftChannel, rightChannel);

    case AudioAnalysisMode::ONSET_DETECTION:
        return analyzeOnsets(leftChannel, rightChannel);

    default:
        return analyzeRMS(leftChannel, rightChannel);
    }
}

AudioFeatures AudioAnalyzer::analyzeRMS(const std::vector<float>& left, const std::vector<float>& right) {
    AudioFeatures features;

    // Calculate raw RMS
    float rawLeftRMS = calculateRMS(left);
    float rawRightRMS = calculateRMS(right);

    if (enableAutoScaling) {
        // Update maximum observed RMS (with slow adaptation)
        float currentMaxRMS = std::max(rawLeftRMS, rawRightRMS);
        if (currentMaxRMS > maxObservedRMS) {
            // Quickly adapt upward for louder audio
            maxObservedRMS = maxObservedRMS * (1.0f - adaptationRate * 10) + currentMaxRMS * (adaptationRate * 10);
        }
        else {
            // Slowly adapt downward to prevent getting stuck at high levels
            maxObservedRMS = maxObservedRMS * (1.0f - adaptationRate) + currentMaxRMS * adaptationRate;
        }

        // Ensure minimum threshold to avoid division by zero
        maxObservedRMS = std::max(maxObservedRMS, 0.001f);

        // Normalize to 0-1 range
        features.leftRMS = std::min(rawLeftRMS / maxObservedRMS, 1.0f);
        features.rightRMS = std::min(rawRightRMS / maxObservedRMS, 1.0f);
    }
    else {
        // Use raw values
        features.leftRMS = rawLeftRMS;
        features.rightRMS = rawRightRMS;
    }

    features.overallEnergy = features.leftRMS + features.rightRMS;
    return features;
}

AudioFeatures AudioAnalyzer::analyzeFFT(const std::vector<float>& left, const std::vector<float>& right) {
    AudioFeatures features;

    // Combine stereo to mono for FFT (you could also analyze separately)
    std::vector<float> mono(left.size());
    for (size_t i = 0; i < left.size() && i < right.size(); i++) {
        mono[i] = (left[i] + right[i]) * 0.5f;
    }

    // Perform FFT
    performFFT(mono, features.fftMagnitudes);

    // Generate frequency bins
    features.fftBins.resize(features.fftMagnitudes.size());
    for (size_t i = 0; i < features.fftBins.size(); i++) {
        features.fftBins[i] = (float)i * sampleRate / (2.0f * features.fftMagnitudes.size());
    }

    // Calculate overall energy from FFT
    features.overallEnergy = std::accumulate(features.fftMagnitudes.begin(),
        features.fftMagnitudes.end(), 0.0f) / features.fftMagnitudes.size();

    return features;
}

AudioFeatures AudioAnalyzer::analyzePeaks(const std::vector<float>& left, const std::vector<float>& right) {
    AudioFeatures features;

    // Combine channels
    std::vector<float> combined(left.size());
    for (size_t i = 0; i < left.size() && i < right.size(); i++) {
        combined[i] = (left[i] + right[i]) * 0.5f;
    }

    features.peakDetected = detectPeak(combined, peakThreshold);
    if (features.peakDetected) {
        features.peakMagnitude = *std::max_element(combined.begin(), combined.end());
    }

    features.overallEnergy = calculateRMS(combined);

    return features;
}

AudioFeatures AudioAnalyzer::analyzeSpectralCentroid(const std::vector<float>& left, const std::vector<float>& right) {
    AudioFeatures features;

    // First perform FFT
    std::vector<float> mono(left.size());
    for (size_t i = 0; i < left.size() && i < right.size(); i++) {
        mono[i] = (left[i] + right[i]) * 0.5f;
    }

    performFFT(mono, features.fftMagnitudes);

    // Generate frequency bins
    std::vector<float> frequencies(features.fftMagnitudes.size());
    for (size_t i = 0; i < frequencies.size(); i++) {
        frequencies[i] = (float)i * sampleRate / (2.0f * frequencies.size());
    }

    features.spectralCentroid = calculateSpectralCentroid(features.fftMagnitudes, frequencies);
    features.overallEnergy = calculateRMS(mono);

    return features;
}

AudioFeatures AudioAnalyzer::analyzeOnsets(const std::vector<float>& left, const std::vector<float>& right) {
    AudioFeatures features;

    // Perform FFT first
    std::vector<float> mono(left.size());
    for (size_t i = 0; i < left.size() && i < right.size(); i++) {
        mono[i] = (left[i] + right[i]) * 0.5f;
    }

    std::vector<float> currentSpectrum;
    performFFT(mono, currentSpectrum);

    features.onsetDetected = detectOnset(currentSpectrum, previousSpectrum);

    if (features.onsetDetected) {
        features.onsetStrength = calculateRMS(mono);
    }

    // Store current spectrum for next frame
    previousSpectrum = currentSpectrum;
    features.overallEnergy = calculateRMS(mono);

    return features;
}

// Helper function implementations
float AudioAnalyzer::calculateRMS(const std::vector<float>& samples) {
    if (samples.empty()) return 0.0f;

    float sum = 0.0f;
    for (float sample : samples) {
        sum += sample * sample;
    }
    return sqrt(sum / samples.size());
}

void AudioAnalyzer::performFFT(const std::vector<float>& input, std::vector<float>& magnitudes) {
    // Simple DFT implementation (for production, use FFTW or similar)
    // This is a basic implementation for demonstration

    int N = std::min((int)input.size(), bufferSize);
    magnitudes.resize(N / 2);

    for (int k = 0; k < N / 2; k++) {
        float real = 0.0f, imag = 0.0f;

        for (int n = 0; n < N; n++) {
            float angle = -2.0f * PI * k * n / N;
            real += input[n] * cos(angle);
            imag += input[n] * sin(angle);
        }

        magnitudes[k] = sqrt(real * real + imag * imag) / N;
    }
}

float AudioAnalyzer::calculateSpectralCentroid(const std::vector<float>& magnitudes,
    const std::vector<float>& frequencies) {
    float weightedSum = 0.0f;
    float magnitudeSum = 0.0f;

    for (size_t i = 0; i < magnitudes.size() && i < frequencies.size(); i++) {
        weightedSum += frequencies[i] * magnitudes[i];
        magnitudeSum += magnitudes[i];
    }

    return (magnitudeSum > 0.0f) ? weightedSum / magnitudeSum : 0.0f;
}

bool AudioAnalyzer::detectPeak(const std::vector<float>& samples, float threshold) {
    if (samples.empty()) return false;

    float currentMax = *std::max_element(samples.begin(), samples.end());
    float previousMax = previousSamples.empty() ? 0.0f :
        *std::max_element(previousSamples.begin(), previousSamples.end());

    previousSamples = samples; // Store for next comparison

    return (currentMax > threshold) && (currentMax > previousMax * 1.2f);
}

bool AudioAnalyzer::detectOnset(const std::vector<float>& currentSpectrum,
    const std::vector<float>& previousSpectrum) {
    if (currentSpectrum.size() != previousSpectrum.size() || currentSpectrum.empty()) {
        return false;
    }

    float spectralDifference = 0.0f;
    for (size_t i = 0; i < currentSpectrum.size(); i++) {
        float diff = currentSpectrum[i] - previousSpectrum[i];
        if (diff > 0) { // Only positive differences (spectral flux)
            spectralDifference += diff;
        }
    }

    return spectralDifference > onsetThreshold;
}

std::vector<std::string> AudioAnalyzer::getModeNames() const {
    return {
        "RMS Amplitude",
        "FFT Spectrum",
        "Peak Detection",
        "Spectral Centroid",
        "Onset Detection"
    };
}

void AudioAnalyzer::nextMode() {
    int current = static_cast<int>(currentMode);
    current = (current + 1) % 5; // Number of modes
    setMode(static_cast<AudioAnalysisMode>(current));
}

void AudioAnalyzer::previousMode() {
    int current = static_cast<int>(currentMode);
    current = (current - 1 + 5) % 5; // Number of modes
    setMode(static_cast<AudioAnalysisMode>(current));
}