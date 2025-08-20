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

    // Initialize spectrum data
    smoothedSpectrum.resize(bufferSize / 2, 0.0f);
    spectrumBars.resize(bufferSize / 2, 0.0f);
    circularSpectrum.resize(bufferSize / 2);
    spectrumHistory.reserve(HISTORY_SIZE);

    // Initialize particles
    particles.reserve(100);
    particleLifetimes.reserve(100);
    particleColors.reserve(100);

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

    if (currentMode == VisualizationMode::WAVE_RMS) {
        updateWaveformData(leftChannel, rightChannel);
    }

    updateRMSValues(features);
    updateSpectrumData(features);
    updatePeakData(features);
    updateOnsetData(features);
    updateCentroidData(features);
    updateParticles(currentOnsetDetected);
    updateVisualState();
}

// Clean draw method - no parameters needed, all data is pre-processed
void AudioVisualizer::draw() {
    ofPushMatrix();
    ofTranslate(visualParams.position);
    ofScale(visualParams.scale);

    // Set default colors
    ofSetColor(visualParams.primaryColor);

    switch (currentMode) {
    case VisualizationMode::WAVE_RMS:
        drawWaveformRMS();
        break;
    case VisualizationMode::SPECTRUM_BARS:
        drawSpectrumBars();
        break;
    case VisualizationMode::CIRCULAR_SPECTRUM:
        drawCircularSpectrum();
        break;
    case VisualizationMode::PEAK_PULSES:
        drawPeakPulses();
        break;
    case VisualizationMode::CENTROID_WAVE:
        drawCentroidWave();
        break;
    case VisualizationMode::ONSET_PARTICLES:
        drawOnsetParticles();
        break;

    case VisualizationMode::COMBINED_VIEW:
        drawCombinedView();
        break;
    }

    // Draw labels if enabled
    if (visualParams.showLabels) {
        std::string modeLabel = getModeNames()[static_cast<int>(currentMode)];
        drawLabel(modeLabel, glm::vec2(10, visualParams.size.y - 20));
    }

    ofPopMatrix();
}

// ============================================================================
// UPDATE METHODS - All data processing happens here
// ============================================================================

void AudioVisualizer::updateRMSValues(const AudioFeatures& features) {
    float factor = visualParams.smoothing;

    smoothedRMS = smoothedRMS * (1.0f - factor) + features.leftRMS * factor;
    leftRMS = leftRMS * (1.0f - factor) + features.leftRMS * factor;
    rightRMS = rightRMS * (1.0f - factor) + features.rightRMS * factor;
}

void AudioVisualizer::updateSpectrumData(const AudioFeatures& features) {
    if (features.fftMagnitudes.empty()) return;

    float factor = visualParams.smoothing;

    // Update smoothed spectrum
    smoothedSpectrum.resize(features.fftMagnitudes.size());
    for (size_t i = 0; i < features.fftMagnitudes.size(); ++i) {
        smoothedSpectrum[i] = smoothedSpectrum[i] * (1.0f - factor) + features.fftMagnitudes[i] * factor;
    }

    // Process spectrum bars (ready-to-render bar heights)
    spectrumBars.resize(features.fftMagnitudes.size());
    float maxHeight = visualParams.size.y * 0.8f;
    for (size_t i = 0; i < features.fftMagnitudes.size(); ++i) {
        spectrumBars[i] = smoothedSpectrum[i] * maxHeight * visualParams.sensitivity;
    }

    // Process circular spectrum (convert to polar coordinates)
    circularSpectrum.resize(features.fftMagnitudes.size());
    float baseRadius = std::min(visualParams.size.x, visualParams.size.y) * 0.2f;
    float maxRadius = std::min(visualParams.size.x, visualParams.size.y) * 0.4f;

    for (size_t i = 0; i < features.fftMagnitudes.size(); ++i) {
        float angle = TWO_PI * i / features.fftMagnitudes.size();
        float magnitude = smoothedSpectrum[i] * visualParams.sensitivity;
        float radius = baseRadius + magnitude * (maxRadius - baseRadius);

        circularSpectrum[i] = polarToCartesian(radius, angle);
    }
}

void AudioVisualizer::updatePeakData(const AudioFeatures& features) {
    currentPeakDetected = features.peakDetected;

    // Update pulse intensity
    if (currentPeakDetected) {
        pulseIntensity = std::min(1.0f, features.peakMagnitude);
    }
    else {
        pulseIntensity *= 0.95f; // Decay
    }
}

void AudioVisualizer::updateOnsetData(const AudioFeatures& features) {
    currentOnsetDetected = features.onsetDetected;
}

void AudioVisualizer::updateCentroidData(const AudioFeatures& features) {
    float factor = visualParams.smoothing;
    smoothedCentroid = smoothedCentroid * (1.0f - factor) + features.spectralCentroid * factor;

    // Normalize centroid for color mapping (0-1 range)
    normalizedCentroid = smoothedCentroid / (sampleRate * 0.5f);
    normalizedCentroid = ofClamp(normalizedCentroid, 0.0f, 1.0f);

    // Pre-calculate color for centroid visualization
    centroidColor = getColorFromFrequency(normalizedCentroid);
}

void AudioVisualizer::updateParticles(bool addNew) {
    // Update existing particles
    for (size_t i = 0; i < particleLifetimes.size(); ++i) {
        particleLifetimes[i] -= 0.02f;
        if (particleLifetimes[i] <= 0) {
            particles.erase(particles.begin() + i);
            particleLifetimes.erase(particleLifetimes.begin() + i);
            particleColors.erase(particleColors.begin() + i);
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
        particleColors.push_back(visualParams.primaryColor);
    }
}

void AudioVisualizer::updateVisualState() {
    // Update rotation for circular visualizations
    rotationAngle += smoothedRMS * 2.0f;
    if (rotationAngle > TWO_PI) rotationAngle -= TWO_PI;

    // Update current color based on audio characteristics
    currentColor = getColorFromAmplitude(smoothedRMS);
}

void AudioVisualizer::updateWaveformData(const std::vector<float>& leftChannel,
    const std::vector<float>& rightChannel) {
    // Downsample waveform data for efficient rendering
    int targetSize = std::min((int)visualParams.size.x, 500); // Limit points for performance

    if (!leftChannel.empty()) {
        processedLeftChannel = downsampleForVisualization(leftChannel, targetSize);
    }

    if (!rightChannel.empty()) {
        processedRightChannel = downsampleForVisualization(rightChannel, targetSize);
    }

    // Update waves if we have audio data
    if (!leftChannel.empty() && waves.size() >= 2) {
        waves[0].setWaveform(leftChannel);   // Left
        if (!rightChannel.empty()) {
            waves[1].setWaveform(rightChannel);  // Right
        }
    }
}

// ============================================================================
// DRAW METHODS - Only rendering, no data processing
// ============================================================================

void AudioVisualizer::drawSpectrumBars() {
    if (spectrumBars.empty()) return;

    ofPushStyle();

    float barWidth = visualParams.size.x / spectrumBars.size();

    for (size_t i = 0; i < spectrumBars.size(); ++i) {
        float height = spectrumBars[i];

        // Color based on frequency
        float freq = (float)i / spectrumBars.size();
        ofColor barColor = getColorFromFrequency(freq);
        ofSetColor(barColor);

        ofDrawRectangle(i * barWidth, visualParams.size.y - height, barWidth - 1, height);
    }

    ofPopStyle();
}

void AudioVisualizer::drawCircularSpectrum() {
    if (circularSpectrum.empty()) return;

    ofPushStyle();

    glm::vec2 center = visualParams.size * 0.5f;

    ofPushMatrix();
    ofTranslate(center);
    ofRotate(ofRadToDeg(rotationAngle));

    for (size_t i = 0; i < circularSpectrum.size(); ++i) {
        glm::vec2 point = circularSpectrum[i];
        float magnitude = smoothedSpectrum[i] * visualParams.sensitivity;

        // Color based on magnitude
        ofColor pointColor = getColorFromAmplitude(magnitude);
        ofSetColor(pointColor);

        ofDrawCircle(point.x, point.y, 2);

        // Draw line from center
        ofSetColor(pointColor, 100);
        ofDrawLine(0, 0, point.x, point.y);
    }

    ofPopMatrix();
    ofPopStyle();
}

void AudioVisualizer::drawPeakPulses() {
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

void AudioVisualizer::drawCentroidWave() {
    if (waves.empty()) return;

    ofPushStyle();

    // Use pre-calculated centroid color
    ofSetColor(centroidColor);

    ofPushMatrix();
    ofTranslate(0, visualParams.size.y * 0.5f);
    waves[0].setAmplitude(smoothedRMS * visualParams.sensitivity * 100);
    waves[0].draw();
    ofPopMatrix();

    ofPopStyle();
}

void AudioVisualizer::drawOnsetParticles() {
    ofPushStyle();

    // Draw existing particles using pre-processed data
    for (size_t i = 0; i < particles.size(); ++i) {
        float life = particleLifetimes[i];
        if (life > 0) {
            ofColor particleColor = particleColors[i];
            particleColor.a = 255 * life;
            ofSetColor(particleColor);
            ofDrawCircle(particles[i], 5 * life);
        }
    }

    ofPopStyle();
}

void AudioVisualizer::drawWaveformRMS() {
    if (processedLeftChannel.empty()) return;

    ofPushStyle();
    ofSetColor(visualParams.primaryColor);
    ofNoFill();
    ofSetLineWidth(2);

    // Draw left channel using processed data
    ofBeginShape();
    for (size_t i = 0; i < processedLeftChannel.size(); ++i) {
        float x = (float)i / processedLeftChannel.size() * visualParams.size.x;
        float y = visualParams.size.y * 0.25f + processedLeftChannel[i] * 100 * visualParams.sensitivity;
        ofVertex(x, y);
    }
    ofEndShape();

    // Draw right channel if available
    if (!processedRightChannel.empty()) {
        ofSetColor(visualParams.secondaryColor);
        ofBeginShape();
        for (size_t i = 0; i < processedRightChannel.size(); ++i) {
            float x = (float)i / processedRightChannel.size() * visualParams.size.x;
            float y = visualParams.size.y * 0.75f + processedRightChannel[i] * 100 * visualParams.sensitivity;
            ofVertex(x, y);
        }
        ofEndShape();
    }

    ofPopStyle();
}

void AudioVisualizer::drawCombinedView() {
    // Combine multiple visualizations in one view
    ofPushMatrix();

    // Top half: spectrum
    ofPushMatrix();
    ofScale(1.0f, 0.4f);
    drawSpectrumBars();
    ofPopMatrix();

    // Bottom half: waveform
    ofPushMatrix();
    ofTranslate(0, visualParams.size.y * 0.5f);
    ofScale(1.0f, 0.5f);
    drawWaveformRMS();
    ofPopMatrix();

    ofPopMatrix();
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

std::vector<float> AudioVisualizer::downsampleForVisualization(const std::vector<float>& input, int targetSize) {
    if (input.empty() || targetSize <= 0) return {};
    if (input.size() <= targetSize) return input;

    std::vector<float> output;
    output.reserve(targetSize);

    float ratio = (float)input.size() / targetSize;
    for (int i = 0; i < targetSize; ++i) {
        int index = (int)(i * ratio);
        if (index < input.size()) {
            output.push_back(input[index]);
        }
    }

    return output;
}

void AudioVisualizer::smoothValue(float& target, float newValue, float smoothing) {
    target = target * (1.0f - smoothing) + newValue * smoothing;
}

glm::vec2 AudioVisualizer::polarToCartesian(float radius, float angle) {
    return glm::vec2(cos(angle) * radius, sin(angle) * radius);
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
    waves.resize(3);  // Left, right, and combined
    for (auto& wave : waves) {
        wave.setup();
    }
}

void AudioVisualizer::previousMode() {
    int currentInt = static_cast<int>(currentMode);
    int maxMode = static_cast<int>(VisualizationMode::COMBINED_VIEW);
    currentInt = (currentInt - 1 + maxMode + 1) % (maxMode + 1);
    setMode(static_cast<VisualizationMode>(currentInt));
    waves.resize(3);  // Left, right, and combined
    for (auto& wave : waves) {
        wave.setup();
    }
}

std::vector<std::string> AudioVisualizer::getModeNames() const {
    return {
        "Wave RMS",
        "Spectrum Bars",
        "Circular Spectrum",
        "Peak Pulses",
        "Centroid Wave",
        "Onset Particles",
        "Combined View"
    };
}

void AudioVisualizer::reset() {
    smoothedRMS = 0.0f;
    leftRMS = 0.0f;
    rightRMS = 0.0f;
    smoothedCentroid = 0.0f;
    normalizedCentroid = 0.0f;
    pulseIntensity = 0.0f;
    rotationAngle = 0.0f;
    currentPeakDetected = false;
    currentOnsetDetected = false;

    particles.clear();
    particleLifetimes.clear();
    particleColors.clear();
    processedLeftChannel.clear();
    processedRightChannel.clear();

    std::fill(smoothedSpectrum.begin(), smoothedSpectrum.end(), 0.0f);
    std::fill(spectrumBars.begin(), spectrumBars.end(), 0.0f);
}

void AudioVisualizer::addToHistory(float value, std::vector<float>& history) {
    history.push_back(value);
    if (history.size() > HISTORY_SIZE) {
        history.erase(history.begin());
    }
}
