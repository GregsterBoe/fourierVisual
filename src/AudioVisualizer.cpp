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

void AudioVisualizer::update(const AudioFeatures& features) {
    currentFeatures = features;

    updateParticles(features.onsetDetected);
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

    ofPopMatrix();
}


// ============================================================================
// DRAW METHODS - Only rendering, no data processing
// ============================================================================

void AudioVisualizer::drawWaveformRMS() {
    if (currentFeatures.leftChannelViz.empty()) return;

    ofPushStyle();
    ofSetColor(visualParams.primaryColor);
    ofNoFill();
    ofSetLineWidth(2);

    // Draw left channel using pre-processed visualization data
    ofBeginShape();
    for (size_t i = 0; i < currentFeatures.leftChannelViz.size(); ++i) {
        float x = (float)i / currentFeatures.leftChannelViz.size() * visualParams.size.x;
        float y = visualParams.size.y * 0.25f + currentFeatures.leftChannelViz[i] * 100 * visualParams.sensitivity;
        ofVertex(x, y);
    }
    ofEndShape();

    // Draw right channel if available
    if (!currentFeatures.rightChannelViz.empty()) {
        ofSetColor(visualParams.secondaryColor);
        ofBeginShape();
        for (size_t i = 0; i < currentFeatures.rightChannelViz.size(); ++i) {
            float x = (float)i / currentFeatures.rightChannelViz.size() * visualParams.size.x;
            float y = visualParams.size.y * 0.75f + currentFeatures.rightChannelViz[i] * 100 * visualParams.sensitivity;
            ofVertex(x, y);
        }
        ofEndShape();
    }

    ofPopStyle();
}

void AudioVisualizer::drawSpectrumBars() {
    if (currentFeatures.fftMagnitudes.empty()) return;

    ofPushStyle();
    auto magnitudes = currentFeatures.logFrequencyBands;
    std::reverse(magnitudes.begin(), magnitudes.end());

    float initialBarWidth = visualParams.size.x / magnitudes.size();
    int groupSize = 1;

    // Calculate how many bins to group together to meet minimum width
    while (initialBarWidth * groupSize < MIN_BAR_WIDTH && groupSize < magnitudes.size()) {
        groupSize *= 2;
    }

    // Calculate final parameters
    int numBars = magnitudes.size() / groupSize;
    float finalBarWidth = visualParams.size.x / numBars;

    // Calculate dominant frequency position for visualization
    float dominantFreqPosition = getDominantFrequencyPosition(currentFeatures.dominantFrequency, numBars);

    for (int i = 0; i < numBars; ++i) {
        // Accumulate neighboring heights
        float accumulatedHeight = 0.0f;
        int actualGroupSize = 0;
        for (int j = 0; j < groupSize && (i * groupSize + j) < magnitudes.size(); ++j) {
            accumulatedHeight += magnitudes[i * groupSize + j];
            actualGroupSize++;
        }

        // Average the accumulated values
        float height = (accumulatedHeight / actualGroupSize) * visualParams.sensitivity * visualParams.size.y;

        // Clean position-based color transition (0.0 = left, 1.0 = right)
        float position = (float)i / (numBars - 1);
        ofColor barColor = getSpectrumColor(position);

        // ENHANCEMENT 1: Highlight bars near dominant frequency
        if (currentFeatures.melodyConfidence > 0.3f) { // Only if we're confident
            float distanceFromDominant = abs(i - dominantFreqPosition);
            if (distanceFromDominant <= 2.0f) { // Within 2 bars of dominant frequency
                float proximityFactor = 1.0f - (distanceFromDominant / 2.0f);
                // Boost brightness and add golden tint for dominant frequency area
                barColor = barColor.lerp(ofColor::gold, proximityFactor * 0.6f);
                barColor.setBrightness(barColor.getBrightness() * (1.0f + proximityFactor * 0.5f));

                // ENHANCEMENT 2: Make dominant frequency bars taller
                if (distanceFromDominant <= 1.0f) {
                    height *= (1.0f + proximityFactor * currentFeatures.melodyConfidence * 0.3f);
                }
            }
        }

        ofSetColor(barColor);
        float x = visualParams.position.x + i * finalBarWidth;
        float y = visualParams.position.y + visualParams.size.y - height;
        ofDrawRectangle(x, y, finalBarWidth - 1, height);
    }

    ofPopStyle();
}

// Helper function to convert frequency to bar position
float AudioVisualizer::getDominantFrequencyPosition(float frequency, int numBars) {
    if (frequency <= 0) return -1;

    // Assuming your frequency bands are logarithmically spaced
    // You'll need to adapt this to match your actual frequency mapping
    float minFreq = 20.0f;   // Lowest frequency in your spectrum
    float maxFreq = 20000.0f; // Highest frequency in your spectrum

    // Logarithmic mapping to match your log frequency bands
    float logFreq = log(frequency);
    float logMin = log(minFreq);
    float logMax = log(maxFreq);

    float normalizedPosition = (logFreq - logMin) / (logMax - logMin);

    // Since you reverse the magnitudes, we need to reverse the position too
    normalizedPosition = 1.0f - normalizedPosition;

    return normalizedPosition * (numBars - 1);
}


void AudioVisualizer::drawCircularSpectrum() {
    if (currentFeatures.circularSpectrum.empty()) return;

    ofPushStyle();

    glm::vec2 center = visualParams.size * 0.5f;

    ofPushMatrix();
    ofTranslate(center);
    ofRotate(ofRadToDeg(currentFeatures.rotationAngle));

    for (size_t i = 0; i < currentFeatures.circularSpectrum.size(); ++i) {
        glm::vec2 point = currentFeatures.circularSpectrum[i];

        // Use magnitude from smoothed spectrum for color
        float magnitude = 0.0f;
        if (i < currentFeatures.smoothedSpectrum.size()) {
            magnitude = currentFeatures.smoothedSpectrum[i] * visualParams.sensitivity;
        }

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
    float pulseRadius = baseRadius + currentFeatures.pulseIntensity * 100 * visualParams.sensitivity;

    // Main pulse circle
    ofColor pulseColor = visualParams.primaryColor;
    pulseColor.a = 255 * (1.0f - currentFeatures.pulseIntensity);
    ofSetColor(pulseColor);
    ofNoFill();
    ofSetLineWidth(3);
    ofDrawCircle(center, pulseRadius);

    // Inner stable circle
    ofSetColor(visualParams.primaryColor);
    ofFill();
    ofDrawCircle(center, baseRadius * currentFeatures.smoothedRMS * visualParams.sensitivity);

    ofPopStyle();
}

void AudioVisualizer::drawCentroidWave() {
    if (waves.empty()) return;

    ofPushStyle();

    // Use pre-calculated centroid color from features
    ofSetColor(currentFeatures.centroidColor);

    ofPushMatrix();
    ofTranslate(0, visualParams.size.y * 0.5f);
    waves[0].setAmplitude(currentFeatures.smoothedRMS * visualParams.sensitivity * 100);
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


void AudioVisualizer::updateParticles(bool addNew) {
    // Update existing particles
    for (size_t i = 0; i < particleLifetimes.size(); ++i) {
        particleLifetimes[i] -= 0.02f;  // Decay rate

        // Remove dead particles
        if (particleLifetimes[i] <= 0) {
            particles.erase(particles.begin() + i);
            particleLifetimes.erase(particleLifetimes.begin() + i);
            particleColors.erase(particleColors.begin() + i);
            --i;  // Adjust index after removal
        }
    }

    // Add new particle on onset
    if (addNew && particles.size() < 50) {  // Limit max particles for performance
        glm::vec2 newParticle = {
            ofRandom(visualParams.size.x),
            ofRandom(visualParams.size.y)
        };

        particles.push_back(newParticle);
        particleLifetimes.push_back(1.0f);  // Start with full life
        particleColors.push_back(visualParams.primaryColor);
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

ofColor AudioVisualizer::getSpectrumColor(float position) {
    // Clamp position to [0,1] range
    position = std::max(0.0f, std::min(1.0f, position));

    // Option 1: Classic rainbow spectrum (ROYGBIV)
    if (position < 0.17f) {
        // Red to Orange
        float t = position / 0.17f;
        return ofColor(255, t * 165, 0);
    }
    else if (position < 0.33f) {
        // Orange to Yellow  
        float t = (position - 0.17f) / 0.16f;
        return ofColor(255, 165 + t * 90, 0);
    }
    else if (position < 0.5f) {
        // Yellow to Green
        float t = (position - 0.33f) / 0.17f;
        return ofColor(255 - t * 255, 255, 0);
    }
    else if (position < 0.67f) {
        // Green to Cyan
        float t = (position - 0.5f) / 0.17f;
        return ofColor(0, 255, t * 255);
    }
    else if (position < 0.83f) {
        // Cyan to Blue
        float t = (position - 0.67f) / 0.16f;
        return ofColor(0, 255 - t * 255, 255);
    }
    else {
        // Blue to Violet
        float t = (position - 0.83f) / 0.17f;
        return ofColor(t * 128, 0, 255);
    }
}

void AudioVisualizer::reset() {
    // Clear all particle data
    particles.clear();
    particleLifetimes.clear();
    particleColors.clear();
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
