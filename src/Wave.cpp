#include "Wave.h"
#include <algorithm>
#include <numeric>
#include <cmath>

Wave::Wave() {
    waveData = generateSineWave(440.0f, 50.0f, 0.0f, 44100, 1024);
}

Wave::Wave(float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    waveData = generateSineWave(frequency, amplitude, phase, sampleRate, numSamples);
}

Wave::Wave(const std::vector<float>& waveformData) {
    waveData = waveformData;
    invalidateCache();
}

void Wave::setup() {
    waveData = generateSineWave(440.0f, 20.0f, 0.0f, 44100, 1024);
    invalidateCache();
}

void Wave::setup(float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    waveData = generateSineWave(frequency, amplitude, phase, sampleRate, numSamples);
    invalidateCache();
}

// Static waveform generators
std::vector<float> Wave::generateSineWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    std::vector<float> wave(numSamples);
    float dt = 1.0f / sampleRate;
    for (int i = 0; i < numSamples; i++) {
        float t = i * dt;
        wave[i] = amplitude * sin(2 * PI * frequency * t + phase);
    }
    return wave;
}

std::vector<float> Wave::generateSquareWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    std::vector<float> wave(numSamples);
    float dt = 1.0f / sampleRate;
    for (int i = 0; i < numSamples; i++) {
        float t = i * dt;
        float sineVal = sin(2 * PI * frequency * t + phase);
        wave[i] = amplitude * (sineVal >= 0 ? 1.0f : -1.0f);
    }
    return wave;
}

std::vector<float> Wave::generateTriangleWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    std::vector<float> wave(numSamples);
    float dt = 1.0f / sampleRate;
    for (int i = 0; i < numSamples; i++) {
        float t = i * dt;
        float normalizedTime = fmod(frequency * t + phase / (2 * PI), 1.0f);
        if (normalizedTime < 0) normalizedTime += 1.0f;

        if (normalizedTime < 0.5f) {
            wave[i] = amplitude * (4.0f * normalizedTime - 1.0f);
        }
        else {
            wave[i] = amplitude * (3.0f - 4.0f * normalizedTime);
        }
    }
    return wave;
}

std::vector<float> Wave::generateSawtoothWave(float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    std::vector<float> wave(numSamples);
    float dt = 1.0f / sampleRate;
    for (int i = 0; i < numSamples; i++) {
        float t = i * dt;
        float normalizedTime = fmod(frequency * t + phase / (2 * PI), 1.0f);
        if (normalizedTime < 0) normalizedTime += 1.0f;
        wave[i] = amplitude * (2.0f * normalizedTime - 1.0f);
    }
    return wave;
}

std::vector<float> Wave::generateNoiseWave(float amplitude, int numSamples) {
    std::vector<float> wave(numSamples);
    for (int i = 0; i < numSamples; i++) {
        wave[i] = amplitude * (2.0f * ofRandom(1.0f) - 1.0f);
    }
    return wave;
}

// Waveform manipulation
void Wave::setWaveform(const std::vector<float>& newWave) {
    waveData = newWave;
    invalidateCache();
}

void Wave::setWaveform(WaveType type, float frequency, float amplitude, float phase, int sampleRate, int numSamples) {
    switch (type) {
    case WaveType::SINE:
        waveData = generateSineWave(frequency, amplitude, phase, sampleRate, numSamples);
        break;
    case WaveType::SQUARE:
        waveData = generateSquareWave(frequency, amplitude, phase, sampleRate, numSamples);
        break;
    case WaveType::TRIANGLE:
        waveData = generateTriangleWave(frequency, amplitude, phase, sampleRate, numSamples);
        break;
    case WaveType::SAWTOOTH:
        waveData = generateSawtoothWave(frequency, amplitude, phase, sampleRate, numSamples);
        break;
    case WaveType::NOISE:
        waveData = generateNoiseWave(amplitude, numSamples);
        break;
    default:
        waveData = generateSineWave(frequency, amplitude, phase, sampleRate, numSamples);
        break;
    }
    invalidateCache();
}

void Wave::setAmplitude(float newAmplitude) {
    if (waveData.empty()) return;

    float currentMax = getPeak();
    if (currentMax > 0) {
        float scaleFactor = newAmplitude / currentMax;
        for (float& sample : waveData) {
            sample *= scaleFactor;
        }
        invalidateCache();
    }
}

void Wave::normalize() {
    if (waveData.empty()) return;

    float maxVal = *std::max_element(waveData.begin(), waveData.end());
    float minVal = *std::min_element(waveData.begin(), waveData.end());
    float range = std::max(std::abs(maxVal), std::abs(minVal));

    if (range > 0) {
        for (float& sample : waveData) {
            sample /= range;
        }
        invalidateCache();
    }
}

void Wave::scale(float factor) {
    for (float& sample : waveData) {
        sample *= factor;
    }
    invalidateCache();
}

void Wave::smooth(float smoothingFactor) {
    if (waveData.size() < 2) return;

    std::vector<float> smoothed = waveData;
    for (size_t i = 1; i < smoothed.size() - 1; i++) {
        smoothed[i] = waveData[i] * (1.0f - smoothingFactor) +
            (waveData[i - 1] + waveData[i + 1]) * 0.5f * smoothingFactor;
    }
    waveData = smoothed;
    invalidateCache();
}

// Drawing methods
void Wave::draw() {
    draw(0, 0, 1.0f, 1.0f);
}

void Wave::draw(float startX, float startY, float scaleX, float scaleY) {
    if (waveData.empty()) return;

    ofPushStyle();
    ofSetColor(drawParams.color, drawParams.opacity);
    ofSetLineWidth(drawParams.lineWidth);

    if (drawParams.showZeroLine) {
        drawZeroLine(glm::vec2(startX, startY),
            glm::vec2(startX + waveData.size() * scaleX, startY));
    }

    if (drawParams.filled) {
        ofBeginShape();
        ofVertex(startX, startY); // Start at zero line
        for (size_t i = 0; i < waveData.size(); i++) {
            float x = startX + i * scaleX;
            float y = startY + waveData[i] * scaleY;
            ofVertex(x, y);
        }
        ofVertex(startX + waveData.size() * scaleX, startY); // End at zero line
        ofEndShape(true);
    }
    else {
        ofPolyline line;
        for (size_t i = 0; i < waveData.size(); i++) {
            float x = startX + i * scaleX;
            float y = startY + waveData[i] * scaleY;
            line.addVertex(x, y);
        }
        line.draw();
    }

    ofPopStyle();
}

void Wave::draw(const glm::vec2& position, const glm::vec2& size) {
    if (waveData.empty()) return;

    float scaleX = size.x / waveData.size();
    float scaleY = size.y / 2.0f; // Divide by 2 assuming wave is normalized to [-1, 1]
    float startY = position.y + size.y / 2.0f; // Center vertically

    draw(position.x, startY, scaleX, scaleY);
}

void Wave::drawWrapped(float centerX, float centerY, float radius) {
    drawWrapped(glm::vec2(centerX, centerY), radius, 1.0f);
}

void Wave::drawWrapped(const glm::vec2& center, float radius, float radiusModulation) {
    if (waveData.empty()) return;

    ofPushStyle();
    ofSetColor(drawParams.color, drawParams.opacity);
    ofSetLineWidth(drawParams.lineWidth);

    ofPolyline line;
    for (size_t i = 0; i < waveData.size(); i++) {
        float angle = ofMap(i, 0, waveData.size(), 0, TWO_PI);
        float r = radius + waveData[i] * radiusModulation;
        float x = center.x + r * cos(angle);
        float y = center.y + r * sin(angle);
        line.addVertex(x, y);
    }

    if (drawParams.filled) {
        line.close();
        line.draw();
        ofFill();
        line.draw();
    }
    else {
        line.close();
        line.draw();
    }

    ofPopStyle();
}

void Wave::drawVertical(float centerX, float startY, float width, float height) {
    if (waveData.empty()) return;

    ofPushStyle();
    ofSetColor(drawParams.color, drawParams.opacity);
    ofSetLineWidth(drawParams.lineWidth);

    float barWidth = width / waveData.size();

    for (size_t i = 0; i < waveData.size(); i++) {
        float x = centerX - width * 0.5f + i * barWidth;
        float barHeight = std::abs(waveData[i]) * height;

        if (drawParams.filled) {
            ofDrawRectangle(x, startY - barHeight, barWidth - 1, barHeight);
        }
        else {
            ofDrawLine(x + barWidth * 0.5f, startY, x + barWidth * 0.5f, startY - barHeight);
        }
    }

    ofPopStyle();
}

void Wave::drawSpectrum(float startX, float startY, float width, float height, bool logarithmic) {
    if (waveData.empty()) return;

    ofPushStyle();
    ofSetColor(drawParams.color, drawParams.opacity);
    ofSetLineWidth(drawParams.lineWidth);

    float barWidth = width / waveData.size();

    for (size_t i = 0; i < waveData.size(); i++) {
        float magnitude = std::abs(waveData[i]);
        if (logarithmic && magnitude > 0) {
            magnitude = 20.0f * log10(magnitude + 0.001f); // Add small value to avoid log(0)
            magnitude = ofMap(magnitude, -60.0f, 0.0f, 0.0f, 1.0f, true); // Map dB to 0-1
        }

        float barHeight = magnitude * height;
        float x = startX + i * barWidth;

        if (drawParams.filled) {
            ofDrawRectangle(x, startY, barWidth - 1, -barHeight);
        }
        else {
            ofDrawLine(x + barWidth * 0.5f, startY, x + barWidth * 0.5f, startY - barHeight);
        }
    }

    ofPopStyle();
}

void Wave::drawAsParticles(const glm::vec2& position, const glm::vec2& size, float particleSize) {
    if (waveData.empty()) return;

    ofPushStyle();
    ofSetColor(drawParams.color, drawParams.opacity);

    for (size_t i = 0; i < waveData.size(); i++) {
        float x = position.x + (float)i / waveData.size() * size.x;
        float y = position.y + size.y * 0.5f + waveData[i] * size.y * 0.5f;

        float alpha = ofMap(std::abs(waveData[i]), 0, getMaxValue(), 50, 255, true);
        ofSetColor(drawParams.color, alpha);
        ofDrawCircle(x, y, particleSize);
    }

    ofPopStyle();
}

void Wave::drawWithColor(const glm::vec2& position, const glm::vec2& size,
    const ofColor& startColor, const ofColor& endColor) {
    if (waveData.empty()) return;

    ofPushStyle();
    ofSetLineWidth(drawParams.lineWidth);

    ofPolyline line;
    for (size_t i = 0; i < waveData.size(); i++) {
        float x = position.x + (float)i / waveData.size() * size.x;
        float y = position.y + size.y * 0.5f + waveData[i] * size.y * 0.25f;

        float t = (float)i / (waveData.size() - 1);
        ofColor currentColor = startColor.getLerped(endColor, t);
        ofSetColor(currentColor, drawParams.opacity);

        if (i > 0) {
            float prevX = position.x + (float)(i - 1) / waveData.size() * size.x;
            float prevY = position.y + size.y * 0.5f + waveData[i - 1] * size.y * 0.25f;
            ofDrawLine(prevX, prevY, x, y);
        }
    }

    ofPopStyle();
}

// Analysis and utility methods
float Wave::getMaxValue() const {
    if (!cacheValid) updateCache();
    return cachedPeak;
}

float Wave::getMinValue() const {
    if (waveData.empty()) return 0.0f;
    return *std::min_element(waveData.begin(), waveData.end());
}

float Wave::getRMS() const {
    if (!cacheValid) updateCache();
    return cachedRMS;
}

float Wave::getPeak() const {
    if (!cacheValid) updateCache();
    return cachedPeak;
}

float Wave::getValueAt(float position) const {
    if (waveData.empty()) return 0.0f;

    position = ofClamp(position, 0.0f, 1.0f);
    float index = position * (waveData.size() - 1);
    return interpolateValue(index);
}

std::vector<float> Wave::getDownsampled(size_t targetSize) const {
    if (waveData.empty() || targetSize == 0) return {};
    if (targetSize >= waveData.size()) return waveData;

    std::vector<float> downsampled(targetSize);
    for (size_t i = 0; i < targetSize; i++) {
        float position = (float)i / (targetSize - 1);
        downsampled[i] = getValueAt(position);
    }
    return downsampled;
}

std::vector<float> Wave::getSmoothed(float factor) const {
    std::vector<float> smoothed = waveData;
    Wave tempWave(smoothed);
    tempWave.smooth(factor);
    return tempWave.getWaveform();
}

// Private helper methods
void Wave::invalidateCache() {
    cacheValid = false;
}

void Wave::updateCache() const {
    if (waveData.empty()) {
        cachedRMS = 0.0f;
        cachedPeak = 0.0f;
        cacheValid = true;
        return;
    }

    // Calculate RMS
    float sum = 0.0f;
    float maxAbs = 0.0f;

    for (float sample : waveData) {
        sum += sample * sample;
        maxAbs = std::max(maxAbs, std::abs(sample));
    }

    cachedRMS = sqrt(sum / waveData.size());
    cachedPeak = maxAbs;
    cacheValid = true;
}

float Wave::interpolateValue(float position) const {
    if (waveData.empty()) return 0.0f;

    int index1 = (int)floor(position);
    int index2 = (int)ceil(position);

    index1 = ofClamp(index1, 0, (int)waveData.size() - 1);
    index2 = ofClamp(index2, 0, (int)waveData.size() - 1);

    if (index1 == index2) {
        return waveData[index1];
    }

    float fraction = position - index1;
    return ofLerp(waveData[index1], waveData[index2], fraction);
}

void Wave::drawZeroLine(const glm::vec2& start, const glm::vec2& end) const {
    ofPushStyle();
    ofSetColor(drawParams.zeroLineColor);
    ofSetLineWidth(1);
    ofDrawLine(start, end);
    ofPopStyle();
}