#include "AudioVisualizer.h"
#include <algorithm>
#include <cmath>

AudioVisualizer::AudioVisualizer()
    : currentMode(VisualizationMode::WAVE_RMS)
    , sampleRate(44100)
    , bufferSize(512)
{
}

AudioVisualizer::~AudioVisualizer() {
}

void AudioVisualizer::setup(int sr, int bs) {
    sampleRate = sr;
    bufferSize = bs;

    // Initialize waves for different visualizations
    waves.resize(3);  // Left, right, and combined
    for (auto& wave : waves) {
        wave.setup();
    }

    // Initialize spectrum history
    smoothedSpectrum.resize(bufferSize / 2, 0.0f);
    spectrumHistory.reserve(HISTORY_SIZE);

    // Initialize particles
    particles.reserve(100);
    particleLifetimes.reserve(100);

    reset();
}

void AudioVisualizer::setMode(VisualizationMode mode) {
    if (currentMode != mode) {
        currentMode = mode;
        reset(); // Reset state when changing modes
    }
}

void AudioVisualizer::update(const AudioFeatures& features,
    const std::vector<float>& leftChannel,
    const std::vector<float>& rightChannel) {
    updateSmoothing(features);
    updateParticles(features.onsetDetected);

    // Update pulse intensity for peak-based visualizations
    if (features.peakDetected) {
        pulseIntensity = std::min(1.0f, features.peakMagnitude);
    }
    else {
        pulseIntensity *= 0.95f; // Decay
    }

    // Update rotation for circular visualizations
    rotationAngle += smoothedRMS * 2.0f;
    if (rotationAngle > TWO_PI) rotationAngle -= TWO_PI;

    // Update waves if we have audio data
    if (!leftChannel.empty() && waves.size() >= 2) {
        waves[0].setWaveform(leftChannel);   // Left
        if (!rightChannel.empty()) {
            waves[1].setWaveform(rightChannel);  // Right
        }
    }
}

void AudioVisualizer::draw(AudioFeatures& features) {
    ofPushMatrix();
    ofTranslate(visualParams.position);
    ofScale(visualParams.scale);

    // Set default colors
    ofSetColor(visualParams.primaryColor);

    // This would typically receive features from the main app
    // For now, we'll show the structure
    AudioFeatures dummyFeatures; // In real use, this comes from update()
    std::vector<float> dummyLeft, dummyRight;



    switch (currentMode) {
    case VisualizationMode::WAVE_RMS:
        drawWaveRMS();
        break;
    case VisualizationMode::SPECTRUM_BARS:
        drawSpectrumBars(dummyFeatures);
        break;
    case VisualizationMode::CIRCULAR_SPECTRUM:
        drawCircularSpectrum(dummyFeatures);
        break;
    case VisualizationMode::PEAK_PULSES:
        drawPeakPulses(dummyFeatures);
        break;
    case VisualizationMode::CENTROID_WAVE:
        drawCentroidWave(dummyFeatures, dummyLeft);
        break;
    case VisualizationMode::ONSET_PARTICLES:
        drawOnsetParticles(dummyFeatures);
        break;
    case VisualizationMode::WAVEFORM_RAW:
        drawRawWaveform(dummyLeft, dummyRight);
        break;
    case VisualizationMode::COMBINED_VIEW:
        drawCombinedView(dummyFeatures, dummyLeft, dummyRight);
        break;
    }

    // Draw labels if enabled
    if (visualParams.showLabels) {
        std::string modeLabel = getModeNames()[static_cast<int>(currentMode)];
        drawLabel(modeLabel, glm::vec2(10, visualParams.size.y - 20));
    }

    ofPopMatrix();
}

void AudioVisualizer::drawWaveRMS() {
    if (waves.empty()) return;

    ofPushStyle();

    // Draw left channel in primary color
    ofSetColor(visualParams.primaryColor, 200);
    ofPushMatrix();
    ofTranslate(0, visualParams.size.y * 0.25f);
    waves[0].setAmplitude(leftRMS * visualParams.sensitivity);
    waves[0].draw(0, 0, defaultScale * 2, defaultScale);
    ofPopMatrix();

    // Draw right channel in secondary color
    if (waves.size() > 1) {
        ofSetColor(visualParams.secondaryColor, 200);
        ofPushMatrix();
        ofTranslate(0, visualParams.size.y * 0.75f);
        waves[1].setAmplitude(rightRMS * visualParams.sensitivity);
        waves[1].draw(0,0, defaultScale * 2, defaultScale);
        ofPopMatrix();
    }

    ofPopStyle();
}

void AudioVisualizer::drawSpectrumBars(const AudioFeatures& features) {
    if (features.fftMagnitudes.empty()) return;

    ofPushStyle();

    float barWidth = visualParams.size.x / features.fftMagnitudes.size();
    float maxHeight = visualParams.size.y * 0.8f;

    for (size_t i = 0; i < features.fftMagnitudes.size() && i < smoothedSpectrum.size(); ++i) {
        float height = smoothedSpectrum[i] * maxHeight * visualParams.sensitivity;

        // Color based on frequency
        float freq = (float)i / features.fftMagnitudes.size();
        ofColor barColor = getColorFromFrequency(freq);
        ofSetColor(barColor);

        ofDrawRectangle(i * barWidth, visualParams.size.y - height, barWidth - 1, height);
    }

    ofPopStyle();
}

void AudioVisualizer::drawCircularSpectrum(const AudioFeatures& features) {
    if (features.fftMagnitudes.empty()) return;

    ofPushStyle();

    glm::vec2 center = visualParams.size * 0.5f;
    float baseRadius = std::min(visualParams.size.x, visualParams.size.y) * 0.2f;
    float maxRadius = std::min(visualParams.size.x, visualParams.size.y) * 0.4f;

    ofPushMatrix();
    ofTranslate(center);
    ofRotate(ofRadToDeg(rotationAngle));

    for (size_t i = 0; i < features.fftMagnitudes.size() && i < smoothedSpectrum.size(); ++i) {
        float angle = TWO_PI * i / features.fftMagnitudes.size();
        float magnitude = smoothedSpectrum[i] * visualParams.sensitivity;
        float radius = baseRadius + magnitude * (maxRadius - baseRadius);

        float x = cos(angle) * radius;
        float y = sin(angle) * radius;

        // Color based on magnitude
        ofColor pointColor = getColorFromAmplitude(magnitude);
        ofSetColor(pointColor);

        ofDrawCircle(x, y, 2);

        // Draw line from center
        ofSetColor(pointColor, 100);
        ofDrawLine(0, 0, x, y);
    }

    ofPopMatrix();
    ofPopStyle();
}

void AudioVisualizer::drawPeakPulses(const AudioFeatures& features) {
    ofPushStyle();

    glm::vec2 center = visualParams.size * 0.5f;
    float baseRadius = 50;
    float pulseRadius = baseRadius + pulseIntensity * 100 * visualParams.sensitivity;

    // Main pulse circle
    ofColor pulseColor = visualParams.primaryColor;
    pulseColor.a = 255 * (1.0f - pulseIntensity);
    ofSetColor(pulseColor);
    ofNoFill();
    ofSetLineWidth(3);
    ofDrawCircle(center, pulseRadius);

    // Inner stable circle
    ofSetColor(visualParams.primaryColor);
    ofFill();
    ofDrawCircle(center, baseRadius * smoothedRMS * visualParams.sensitivity);

    ofPopStyle();
}

void AudioVisualizer::drawCentroidWave(const AudioFeatures& features,
    const std::vector<float>& leftChannel) {
    if (waves.empty() || leftChannel.empty()) return;

    ofPushStyle();

    // Color the wave based on spectral centroid
    ofColor waveColor = getColorFromFrequency(smoothedCentroid / (sampleRate * 0.5f));
    ofSetColor(waveColor);

    ofPushMatrix();
    ofTranslate(0, visualParams.size.y * 0.5f);
    waves[0].setAmplitude(smoothedRMS * visualParams.sensitivity * 100);
    waves[0].draw();
    ofPopMatrix();

    ofPopStyle();
}

void AudioVisualizer::drawOnsetParticles(const AudioFeatures& features) {
    ofPushStyle();

    // Draw existing particles
    for (size_t i = 0; i < particles.size(); ++i) {
        float life = particleLifetimes[i];
        if (life > 0) {
            ofColor particleColor = visualParams.primaryColor;
            particleColor.a = 255 * life;
            ofSetColor(particleColor);
            ofDrawCircle(particles[i], 5 * life);
        }
    }

    ofPopStyle();
}

void AudioVisualizer::drawRawWaveform(const std::vector<float>& leftChannel,
    const std::vector<float>& rightChannel) {
    if (leftChannel.empty()) return;

    ofPushStyle();
    ofSetColor(visualParams.primaryColor);
    ofNoFill();
    ofSetLineWidth(2);

    // Draw left channel
    ofBeginShape();
    for (size_t i = 0; i < leftChannel.size(); ++i) {
        float x = (float)i / leftChannel.size() * visualParams.size.x;
        float y = visualParams.size.y * 0.25f + leftChannel[i] * 100 * visualParams.sensitivity;
        ofVertex(x, y);
    }
    ofEndShape();

    // Draw right channel if available
    if (!rightChannel.empty()) {
        ofSetColor(visualParams.secondaryColor);
        ofBeginShape();
        for (size_t i = 0; i < rightChannel.size(); ++i) {
            float x = (float)i / rightChannel.size() * visualParams.size.x;
            float y = visualParams.size.y * 0.75f + rightChannel[i] * 100 * visualParams.sensitivity;
            ofVertex(x, y);
        }
        ofEndShape();
    }

    ofPopStyle();
}

void AudioVisualizer::drawCombinedView(const AudioFeatures& features,
    const std::vector<float>& leftChannel,
    const std::vector<float>& rightChannel) {
    // Combine multiple visualizations in one view
    ofPushMatrix();

    // Top half: spectrum
    ofPushMatrix();
    ofScale(1.0f, 0.4f);
    drawSpectrumBars(features);
    ofPopMatrix();

    // Bottom half: waveform
    ofPushMatrix();
    ofTranslate(0, visualParams.size.y * 0.5f);
    ofScale(1.0f, 0.5f);
    drawRawWaveform(leftChannel, rightChannel);
    ofPopMatrix();

    ofPopMatrix();
}

// Helper functions implementation
void AudioVisualizer::updateParticles(bool addNew) {
    // Update existing particles
    for (size_t i = 0; i < particleLifetimes.size(); ++i) {
        particleLifetimes[i] -= 0.02f;
        if (particleLifetimes[i] <= 0) {
            particles.erase(particles.begin() + i);
            particleLifetimes.erase(particleLifetimes.begin() + i);
            --i;
        }
    }

    // Add new particle on onset
    if (addNew && particles.size() < 50) {
        glm::vec2 newParticle = {
            ofRandom(visualParams.size.x),
            ofRandom(visualParams.size.y)
        };
        particles.push_back(newParticle);
        particleLifetimes.push_back(1.0f);
    }
}

void AudioVisualizer::updateSmoothing(const AudioFeatures& features) {
    float factor = visualParams.smoothing;

    // Smooth RMS
    smoothedRMS = smoothedRMS * (1.0f - factor) + features.leftRMS * factor;
    leftRMS = leftRMS * (1.0f - factor) + features.leftRMS * factor;
    rightRMS = rightRMS * (1.0f - factor) + features.rightRMS * factor;


    // Smooth spectral centroid
    smoothedCentroid = smoothedCentroid * (1.0f - factor) + features.spectralCentroid * factor;

    // Smooth spectrum
    if (!features.fftMagnitudes.empty()) {
        smoothedSpectrum.resize(features.fftMagnitudes.size());
        for (size_t i = 0; i < features.fftMagnitudes.size(); ++i) {
            smoothedSpectrum[i] = smoothedSpectrum[i] * (1.0f - factor) + features.fftMagnitudes[i] * factor;
        }
    }
}

ofColor AudioVisualizer::getColorFromFrequency(float freq) {
    // Map frequency to color (low = red, mid = green, high = blue)
    freq = ofClamp(freq, 0.0f, 1.0f);
    if (freq < 0.33f) {
        return ofColor(255, freq * 765, 0); // Red to yellow
    }
    else if (freq < 0.67f) {
        return ofColor(255 - (freq - 0.33f) * 765, 255, 0); // Yellow to green
    }
    else {
        return ofColor(0, 255 - (freq - 0.67f) * 765, (freq - 0.67f) * 765); // Green to blue
    }
}

ofColor AudioVisualizer::getColorFromAmplitude(float amplitude) {
    amplitude = ofClamp(amplitude, 0.0f, 1.0f);
    return ofColor(visualParams.primaryColor.r * amplitude,
        visualParams.primaryColor.g * amplitude,
        visualParams.primaryColor.b * amplitude);
}

void AudioVisualizer::drawLabel(const std::string& text, glm::vec2 position) {
    ofPushStyle();
    ofSetColor(255, 200);
    ofDrawBitmapString(text, position);
    ofPopStyle();
}

void AudioVisualizer::nextMode() {
    int currentInt = static_cast<int>(currentMode);
    int maxMode = static_cast<int>(VisualizationMode::COMBINED_VIEW);
    currentInt = (currentInt + 1) % (maxMode + 1);
    setMode(static_cast<VisualizationMode>(currentInt));
}

void AudioVisualizer::previousMode() {
    int currentInt = static_cast<int>(currentMode);
    int maxMode = static_cast<int>(VisualizationMode::COMBINED_VIEW);
    currentInt = (currentInt - 1 + maxMode + 1) % (maxMode + 1);
    setMode(static_cast<VisualizationMode>(currentInt));
}

std::vector<std::string> AudioVisualizer::getModeNames() const {
    return {
        "Wave RMS",
        "Spectrum Bars",
        "Circular Spectrum",
        "Peak Pulses",
        "Centroid Wave",
        "Onset Particles",
        "Raw Waveform",
        "Combined View"
    };
}

void AudioVisualizer::reset() {
    smoothedRMS = 0.0f;
    smoothedCentroid = 0.0f;
    pulseIntensity = 0.0f;
    rotationAngle = 0.0f;
    particles.clear();
    particleLifetimes.clear();
    std::fill(smoothedSpectrum.begin(), smoothedSpectrum.end(), 0.0f);
}

void AudioVisualizer::addToHistory(float value, std::vector<float>& history) {
    history.push_back(value);
    if (history.size() > HISTORY_SIZE) {
        history.erase(history.begin());
    }
}