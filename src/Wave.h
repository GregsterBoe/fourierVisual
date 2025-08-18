#pragma once
#include "ofMain.h"

enum class WaveType {
    SINE,
    COSINE,
    SQUARE,
    TRIANGLE,
    SAWTOOTH,
    NOISE,
    CUSTOM  // For user-provided waveform data
};

class Wave {
public:
    // Constructors
    Wave();
    Wave(float frequency, float amplitude, float phase, int sampleRate, int numSamples);
    Wave(const std::vector<float>& waveformData);  // For real audio data

    // Setup and generation
    void setup();
    void setup(float frequency, float amplitude, float phase, int sampleRate, int numSamples);

    // Waveform generation (static methods for utility)
    static std::vector<float> generateSineWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples);
    static std::vector<float> generateSquareWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples);
    static std::vector<float> generateTriangleWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples);
    static std::vector<float> generateSawtoothWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples);
    static std::vector<float> generateNoiseWave(float amplitude, int numSamples);

    // Waveform manipulation
    void setWaveform(const std::vector<float>& newWave);
    void setWaveform(WaveType type, float frequency, float amplitude, float phase, int sampleRate, int numSamples);
    void setAmplitude(float newAmplitude);
    void setFrequency(float newFrequency, int sampleRate);
    void setPhase(float newPhase, int sampleRate);
    void normalize();
    void scale(float factor);
    void smooth(float smoothingFactor = 0.1f);

    // Drawing methods with more flexibility
    void draw();  // Draw with current parameters
    void draw(float startX, float startY, float scaleX, float scaleY);
    void draw(const glm::vec2& position, const glm::vec2& size);
    void drawWrapped(float centerX, float centerY, float radius);
    void drawWrapped(const glm::vec2& center, float radius, float radiusModulation = 1.0f);
    void drawVertical(float centerX, float startY, float width, float height);
    void drawSpectrum(float startX, float startY, float width, float height, bool logarithmic = false);

    // Advanced drawing modes
    void drawAsParticles(const glm::vec2& position, const glm::vec2& size, float particleSize = 2.0f);
    void drawAs3D(const glm::vec2& position, const glm::vec2& size, float depth = 50.0f);
    void drawWithColor(const glm::vec2& position, const glm::vec2& size, const ofColor& startColor, const ofColor& endColor);

    // Getters
    const std::vector<float>& getWaveform() const { return waveData; }
    std::vector<float>& getWaveform() { return waveData; }
    float getMaxValue() const;
    float getMinValue() const;
    float getRMS() const;
    float getPeak() const;
    size_t getSize() const { return waveData.size(); }
    bool isEmpty() const { return waveData.empty(); }

    // Analysis helpers
    float getValueAt(float position) const;  // position 0-1, interpolated
    std::vector<float> getDownsampled(size_t targetSize) const;
    std::vector<float> getSmoothed(float factor) const;

    // Visualization parameters
    struct DrawParams {
        ofColor color = ofColor::white;
        float lineWidth = 1.0f;
        bool filled = false;
        bool smoothed = true;
        float opacity = 255.0f;
        bool showZeroLine = false;
        ofColor zeroLineColor = ofColor::gray;
    };

    void setDrawParams(const DrawParams& params) { drawParams = params; }
    DrawParams& getDrawParams() { return drawParams; }
    const DrawParams& getDrawParams() const { return drawParams; }

private:
    std::vector<float> waveData;
    DrawParams drawParams;

    // Cached values for performance
    mutable float cachedRMS = -1.0f;
    mutable float cachedPeak = -1.0f;
    mutable bool cacheValid = false;

    // Helper methods
    void invalidateCache();
    void updateCache() const;
    float interpolateValue(float position) const;
    void drawZeroLine(const glm::vec2& start, const glm::vec2& end) const;
};