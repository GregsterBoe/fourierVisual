#include "AudioAnalyzer.h"
#include <algorithm>
#include <numeric>

AudioAnalyzer::AudioAnalyzer()
    : sampleRate(44100)
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



AudioFeatures AudioAnalyzer::analyze(const std::vector<float>&leftChannel,
    const std::vector<float>&rightChannel) {

    AudioFeatures features;

    // Always compute basic RMS (needed by multiple visualizations)
    features.leftRMS = calculateRMS(leftChannel);
    features.rightRMS = calculateRMS(rightChannel);

    // Apply auto-scaling if enabled
    if (enableAutoScaling) {
        float currentMaxRMS = std::max(features.leftRMS, features.rightRMS);
        if (currentMaxRMS > maxObservedRMS) {
            maxObservedRMS = maxObservedRMS * (1.0f - adaptationRate * 10) + currentMaxRMS * (adaptationRate * 10);
        }
        else {
            maxObservedRMS = maxObservedRMS * (1.0f - adaptationRate) + currentMaxRMS * adaptationRate;
        }
        maxObservedRMS = std::max(maxObservedRMS, 0.001f);

        features.leftRMS = std::min(features.leftRMS / maxObservedRMS, 1.0f);
        features.rightRMS = std::min(features.rightRMS / maxObservedRMS, 1.0f);
    }

    // Calculate smoothed RMS for visualizations
    float currentRMS = (features.leftRMS + features.rightRMS) * 0.5f;
    smoothedRMS = smoothedRMS * (1.0f - visualSmoothingFactor) + currentRMS * visualSmoothingFactor;
    features.smoothedRMS = smoothedRMS;

    std::vector<float> mono(leftChannel.size());
    for (size_t i = 0; i < leftChannel.size() && i < rightChannel.size(); i++) {
        mono[i] = (leftChannel[i] + rightChannel[i]) * 0.5f;
    }

    performFFTOptimal(mono, features.fftMagnitudes);

    // Generate frequency bins (same as before)
    features.fftBins.resize(features.fftMagnitudes.size());
    for (size_t i = 0; i < features.fftBins.size(); i++) {
        features.fftBins[i] = (float)i * sampleRate / (2.0f * features.fftMagnitudes.size());
    }

    if (enableMelodyTracking) {
        melodyTracker.updateRange(features.fftMagnitudes, sampleRate);
        features.logFrequencyBands = getLogFrequencyBands(features.fftMagnitudes, bufferSize);
        features.dominantFrequency = melodyTracker.getDominantFrequency();
        features.melodyConfidence = melodyTracker.getMelodyConfidence();
        features.melodyRange = melodyTracker.getCurrentMelodyRange();
    }
    else {
        features.logFrequencyBands = getLogFrequencyBands(features.fftMagnitudes, bufferSize);

    }

    // Update running maximum
    float currentMax = *std::max_element(features.logFrequencyBands.begin(), features.logFrequencyBands.end());
    runningMax = std::max(currentMax, runningMax * maxDecayRate);

    // Normalize using running maximum
    if (runningMax > 0.0f) {
        for (auto& band : features.logFrequencyBands) {
            band = std::min(1.0f, band / runningMax); // Clamp to [0, 1]
        }
    }
    // Calculate energy in different frequency ranges
    calculateFrequencyEnergies(features);

    // Generate frequency bins
    features.fftBins.resize(features.fftMagnitudes.size());
    for (size_t i = 0; i < features.fftBins.size(); i++) {
        features.fftBins[i] = (float)i * sampleRate / (2.0f * features.fftMagnitudes.size());
    }

    // Process smoothed spectrum
    if (smoothedSpectrum.size() != features.fftMagnitudes.size()) {
        smoothedSpectrum.resize(features.fftMagnitudes.size(), 0.0f);
    }

    for (size_t i = 0; i < features.fftMagnitudes.size(); ++i) {
        smoothedSpectrum[i] = smoothedSpectrum[i] * (1.0f - visualSmoothingFactor) +
            features.fftMagnitudes[i] * visualSmoothingFactor;
    }
    features.smoothedSpectrum = smoothedSpectrum;

    // Process circular spectrum (polar coordinates)
    features.circularSpectrum.resize(features.fftMagnitudes.size());
    float baseRadius = std::min(visualizationSize.x, visualizationSize.y) * 0.2f;
    float maxRadius = std::min(visualizationSize.x, visualizationSize.y) * 0.4f;

    for (size_t i = 0; i < features.fftMagnitudes.size(); ++i) {
        float angle = 2.0f * PI * i / features.fftMagnitudes.size();
        float magnitude = smoothedSpectrum[i];
        float radius = baseRadius + magnitude * (maxRadius - baseRadius);

        features.circularSpectrum[i] = glm::vec2(
            cos(angle) * radius,
            sin(angle) * radius
        );
    }

    // Peak detection
    features.peakDetected = detectPeak(mono, peakThreshold);
    if (features.peakDetected) {
        features.peakMagnitude = *std::max_element(mono.begin(), mono.end());
        pulseIntensity = std::min(1.0f, features.peakMagnitude);
    }
    else {
        pulseIntensity *= 0.95f; // Decay
    }
    features.pulseIntensity = pulseIntensity;

    // Spectral centroid
    features.spectralCentroid = calculateSpectralCentroid(features.fftMagnitudes, features.fftBins);
    smoothedCentroid = smoothedCentroid * (1.0f - visualSmoothingFactor) +
        features.spectralCentroid * visualSmoothingFactor;
    features.smoothedCentroid = smoothedCentroid;

    // Normalize centroid for color mapping (0-1 range)
    features.normalizedCentroid = smoothedCentroid / (sampleRate * 0.5f);
    features.normalizedCentroid = std::max(0.0f, std::min(features.normalizedCentroid, 1.0f));

    // Pre-calculate centroid color
    features.centroidColor = getColorFromFrequency(features.normalizedCentroid);

    // Onset detection
    features.onsetDetected = detectOnset(features.fftMagnitudes, previousMagnitudes);
    if (features.onsetDetected) {
        features.onsetStrength = calculateRMS(mono);
    }
    previousMagnitudes = features.fftMagnitudes;

    // Process waveform data for visualization (downsample)
    int targetSize = std::min((int)visualizationSize.x, 500);
    features.leftChannelViz = downsampleForVisualization(leftChannel, targetSize);
    features.rightChannelViz = downsampleForVisualization(rightChannel, targetSize);

    // Update rotation angle for circular visualizations
    rotationAngle += smoothedRMS * 2.0f;
    if (rotationAngle > 2.0f * PI) rotationAngle -= 2.0f * PI;
    features.rotationAngle = rotationAngle;

    // Calculate current color based on audio characteristics
    features.currentColor = getColorFromAmplitude(smoothedRMS);

    // General energy
    features.overallEnergy = smoothedRMS;

    return features;
}

// Add these helper methods to AudioAnalyzer class:

std::vector<float> AudioAnalyzer::downsampleForVisualization(const std::vector<float>& input, int targetSize) {
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

ofColor AudioAnalyzer::getColorFromFrequency(float freq) {
    // Map frequency to color (low = red, mid = green, high = blue)
    freq = std::max(0.0f, std::min(freq, 1.0f));
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

ofColor AudioAnalyzer::getColorFromAmplitude(float amplitude) {
    amplitude = std::max(0.0f, std::min(amplitude, 1.0f));
    // Return a color that scales with amplitude (white base)
    return ofColor(255 * amplitude, 255 * amplitude, 255 * amplitude);
}

float AudioAnalyzer::calculateRMS(const std::vector<float>& samples) {
    if (samples.empty()) return 0.0f;

    float sum = 0.0f;
    for (float sample : samples) {
        sum += sample * sample;
    }
    return sqrt(sum / samples.size());
}

// Fast FFT implementation using Cooley-Tukey algorithm
// Add this as an alternative method to your AudioAnalyzer class

void AudioAnalyzer::performFFTFast(const std::vector<float>& input, std::vector<float>& magnitudes) {
    int N = std::min((int)input.size(), bufferSize);

    // Ensure N is a power of 2 for efficient FFT
    int fftSize = 1;
    while (fftSize < N) fftSize <<= 1;

    // Initialize complex input array
    std::vector<std::complex<float>> fft_input(fftSize, 0.0f);

    // Apply windowing - use Hann window for better frequency resolution
    for (int i = 0; i < N; i++) {
        // Hann window - better than Hamming for frequency analysis
        float window = 0.5f * (1.0f - cos(2.0f * PI * i / (N - 1)));
        fft_input[i] = std::complex<float>(input[i] * window, 0.0f);
    }

    // Perform FFT
    cooleyTukeyFFT(fft_input);

    // Extract magnitudes with minimal processing to preserve selectivity
    int numBins = fftSize / 2;
    magnitudes.resize(numBins);

    for (int i = 0; i < numBins; i++) {
        float magnitude = std::abs(fft_input[i]);

        // Simple, clean normalization
        magnitude = magnitude / fftSize;

        // Skip DC bin (often just noise)
        if (i == 0) {
            magnitude = 0.0f;
        }

        magnitudes[i] = magnitude;
    }

    // Apply minimal, frequency-aware processing
    enhanceFrequencySelectivity(magnitudes);
}

void AudioAnalyzer::enhanceFrequencySelectivity(std::vector<float>& magnitudes) {
    if (magnitudes.empty()) return;

    float sampleRateFloat = (float)sampleRate;
    float binWidth = sampleRateFloat / (magnitudes.size() * 2);

    for (size_t i = 1; i < magnitudes.size(); ++i) {
        float frequency = i * binWidth;
        float& magnitude = magnitudes[i];

        // Apply frequency-specific enhancements
        if (frequency < 60.0f) {
            // Sub-bass: like chill
            magnitude *= 0.5f;
        }
        else if (frequency < 250.0f) {
            // Bass: Slow dooown boi
            magnitude *= 0.4f;
        }
        else if (frequency < 2000.0f) {
            // Mid-range: Slight boost
            magnitude *= 1.1f;
        }
        else if (frequency < 8000.0f) {
            // Upper mids: Increase
            magnitude *= 2.0f;
        }
        else {
            // Highs: puuush
            magnitude *= 4.0f;
        }

        // Apply dynamic range enhancement
        if (magnitude > 0.001f) {
            // Power law to increase dynamic range
            magnitude = pow(magnitude, 0.5f);
        }
    }
}

std::vector<float> AudioAnalyzer::getLogFrequencyBands(const std::vector<float>& magnitudes, int numBands) {
    std::vector<float> logBands(numBands, 0.0f);

    if (magnitudes.empty() || numBands <= 0) return logBands;

    float sampleRateFloat = (float)sampleRate;
    float binWidth = sampleRateFloat / (magnitudes.size() * 2);
    float maxFreq = sampleRateFloat * 0.5f;
    float minFreq = 20.0f;

    // Calculate logarithmic frequency boundaries
    std::vector<float> boundaries(numBands + 1);
    float logMin = log10(minFreq);
    float logMax = log10(maxFreq);
    float logStep = (logMax - logMin) / numBands;

    for (int i = 0; i <= numBands; i++) {
        boundaries[i] = pow(10.0f, logMin + i * logStep);
    }

    // Use different aggregation strategies to preserve separation
    for (int band = 0; band < numBands; band++) {
        float lowFreq = boundaries[band];
        float highFreq = boundaries[band + 1];

        int startBin = std::max(1, (int)(lowFreq / binWidth));
        int endBin = std::min((int)magnitudes.size() - 1, (int)(highFreq / binWidth));

        if (startBin <= endBin) {

            // Strategy 3: Power mean (preserves dynamic range better than arithmetic mean)
            float powerSum = 0.0f;
            int count = 0;
            for (int bin = startBin; bin <= endBin; bin++) {
                powerSum += pow(magnitudes[bin], 2.0f); // Square each value
                count++;
             }
             logBands[band] = count > 0 ? sqrt(powerSum / count) : 0.0f;
        }
    }

    return logBands;
}

// Cooley-Tukey FFT algorithm (add this as a private method)
void AudioAnalyzer::cooleyTukeyFFT(std::vector<std::complex<float>>& data) {
    const int N = data.size();
    if (N <= 1) return;

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }

    // Cooley-Tukey FFT
    for (int len = 2; len <= N; len <<= 1) {
        float ang = 2.0f * PI / len;
        std::complex<float> wlen(cos(ang), sin(ang));

        for (int i = 0; i < N; i += len) {
            std::complex<float> w(1.0f, 0.0f);

            for (int j = 0; j < len / 2; j++) {
                std::complex<float> u = data[i + j];
                std::complex<float> v = data[i + j + len / 2] * w;

                data[i + j] = u + v;
                data[i + j + len / 2] = u - v;

                w *= wlen;
            }
        }
    }
}

// Method 3: Real-valued FFT optimized for audio
// This is more efficient for real input signals
void AudioAnalyzer::performRealFFT(const std::vector<float>& input, std::vector<float>& magnitudes) {
    int N = std::min((int)input.size(), bufferSize);

    // Ensure even size
    if (N % 2 != 0) N--;

    magnitudes.resize(N / 2 + 1); // For real FFT, we get N/2+1 frequency bins

    // Apply windowing
    std::vector<float> windowed(N);
    for (int n = 0; n < N; n++) {
        // Hann window (alternative to Hamming, often better for audio)
        float window = 0.5f * (1.0f - cos(2.0f * PI * n / (N - 1)));
        windowed[n] = input[n] * window;
    }

    // Compute real FFT using DFT
    for (int k = 0; k <= N / 2; k++) {
        float real = 0.0f, imag = 0.0f;

        for (int n = 0; n < N; n++) {
            float angle = -2.0f * PI * k * n / N;
            real += windowed[n] * cos(angle);
            imag += windowed[n] * sin(angle);
        }

        float magnitude = sqrt(real * real + imag * imag);

        // Normalize and scale for visualization
        magnitude = magnitude / N;

        // For DC and Nyquist components, don't double
        if (k > 0 && k < N / 2) {
            magnitude *= 2.0f; // Account for negative frequencies
        }

        // Apply mild logarithmic compression for better visualization
        magnitudes[k] = sqrt(magnitude); // Square root gives good visual balance
    }
}

// Utility method to choose the best FFT method based on buffer size
void AudioAnalyzer::performFFTOptimal(const std::vector<float>& input, std::vector<float>& magnitudes) {
    int N = std::min((int)input.size(), bufferSize);

    // For small buffer sizes, use regular DFT
    if (N < 64) {
        performRealFFT(input, magnitudes);
    }
    // For power-of-2 sizes >= 64, use fast FFT
    else {
        performFFTFast(input, magnitudes);
    }

    applySelectiveSmoothing(magnitudes);
}

void AudioAnalyzer::applySelectiveSmoothing(std::vector<float>& magnitudes) {
    if (magnitudes.size() != smoothedSpectrum.size()) {
        smoothedSpectrum.resize(magnitudes.size(), 0.0f);
    }

    float sampleRateFloat = (float)sampleRate;
    float binWidth = sampleRateFloat / (magnitudes.size() * 2);

    for (size_t i = 0; i < magnitudes.size(); ++i) {
        float frequency = i * binWidth;
        float currentValue = magnitudes[i];
        float& smoothedValue = smoothedSpectrum[i];

        // Adaptive smoothing based on signal characteristics
        float smoothingFactor = 0.1f;

        // Peak-aware smoothing: less smoothing when signal is rising
        if (currentValue > smoothedValue * 1.1f) {
            smoothingFactor *= 0.5f; // Reduce smoothing for rising signals
        }

        smoothedValue = smoothedValue * (1.0f - smoothingFactor) + currentValue * smoothingFactor;
        magnitudes[i] = smoothedValue;
    }
}

void AudioAnalyzer::calculateFrequencyEnergies(AudioFeatures& features) {
    if (features.fftMagnitudes.empty()) return;

    float sampleRateFloat = (float)sampleRate;
    float binWidth = sampleRateFloat / (features.fftMagnitudes.size() * 2);

    float bassSum = 0.0f, midSum = 0.0f, trebleSum = 0.0f;
    int bassCount = 0, midCount = 0, trebleCount = 0;

    for (size_t i = 1; i < features.fftMagnitudes.size(); i++) {
        float frequency = i * binWidth;
        float magnitude = features.fftMagnitudes[i];

        if (frequency < 250.0f) {
            bassSum += magnitude;
            bassCount++;
        }
        else if (frequency < 4000.0f) {
            midSum += magnitude;
            midCount++;
        }
        else {
            trebleSum += magnitude;
            trebleCount++;
        }
    }

    features.bassEnergy = bassCount > 0 ? bassSum / bassCount : 0.0f;
    features.midEnergy = midCount > 0 ? midSum / midCount : 0.0f;
    features.trebleEnergy = trebleCount > 0 ? trebleSum / trebleCount : 0.0f;
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

void MelodyTracker::reset() {
    // Default to common vocal/instrument range
    currentRange = MelodyRange(200.0f, 1500.0f, 440.0f, 0.0f); // A4 as default centroid
    energyHistory.clear();
    previousBands.clear();
}

void MelodyTracker::updateRange(const std::vector<float>& magnitudes, float sampleRate) {
    if (magnitudes.empty()) return;

    float binWidth = sampleRate / (magnitudes.size() * 2);

    // Calculate total energy in melody range for normalization
    float totalMelodyEnergy = calculateMelodyEnergy(magnitudes, binWidth);

    if (totalMelodyEnergy < energyThreshold) {
        // Very quiet, decay confidence but keep range
        currentRange.confidence *= confidenceDecay;
        return;
    }

    // Store energy history
    energyHistory.push_back(totalMelodyEnergy);
    if (energyHistory.size() > historyFrames) {
        energyHistory.erase(energyHistory.begin());
    }

    // Update melody range based on current spectrum
    //updateMelodyRange(magnitudes, binWidth);
}

float MelodyTracker::calculateMelodyCentroid(const std::vector<float>& magnitudes, float binWidth) {
    float weightedSum = 0.0f;
    float magnitudeSum = 0.0f;

    int startBin = std::max(1, (int)(melodyMin / binWidth));
    int endBin = std::min((int)magnitudes.size() - 1, (int)(melodyMax / binWidth));

    // Focus on the most prominent frequencies (peak detection approach)
    for (int i = startBin; i <= endBin; i++) {
        float frequency = i * binWidth;
        float magnitude = magnitudes[i];

        // Emphasize peaks for better note detection
        if (i > 0 && i < magnitudes.size() - 1) {
            // Only consider local maxima and strong signals
            bool isLocalMax = magnitude > magnitudes[i - 1] && magnitude > magnitudes[i + 1];
            if (isLocalMax && magnitude > energyThreshold * 2.0f) {
                // Weight peaks more heavily
                magnitude *= 3.0f;
            }
        }

        weightedSum += frequency * magnitude;
        magnitudeSum += magnitude;
    }

    return (magnitudeSum > 0.0f) ? weightedSum / magnitudeSum : currentRange.centroid;
}

float MelodyTracker::calculateMelodyEnergy(const std::vector<float>& magnitudes, float binWidth) {
    float energy = 0.0f;
    int count = 0;

    int startBin = std::max(1, (int)(melodyMin / binWidth));
    int endBin = std::min((int)magnitudes.size() - 1, (int)(melodyMax / binWidth));

    for (int i = startBin; i <= endBin; i++) {
        energy += magnitudes[i];
        count++;
    }

    return count > 0 ? energy / count : 0.0f;
}

std::pair<float, float> MelodyTracker::findMelodyPeaks(const std::vector<float>& magnitudes, float binWidth) {
    float primaryPeak = 0.0f;
    float secondaryPeak = 0.0f;
    float primaryFreq = currentRange.centroid;
    float secondaryFreq = currentRange.centroid;

    int startBin = std::max(1, (int)(melodyMin / binWidth));
    int endBin = std::min((int)magnitudes.size() - 1, (int)(melodyMax / binWidth));

    // Find the two strongest peaks
    for (int i = startBin + 1; i < endBin - 1; i++) {
        float frequency = i * binWidth;
        float magnitude = magnitudes[i];

        // Check if this is a local maximum
        if (magnitude > magnitudes[i - 1] && magnitude > magnitudes[i + 1] &&
            magnitude > energyThreshold) {

            if (magnitude > primaryPeak) {
                // New primary peak
                secondaryPeak = primaryPeak;
                secondaryFreq = primaryFreq;
                primaryPeak = magnitude;
                primaryFreq = frequency;
            }
            else if (magnitude > secondaryPeak) {
                // New secondary peak
                secondaryPeak = magnitude;
                secondaryFreq = frequency;
            }
        }
    }

    return { primaryFreq, secondaryFreq };
}

void MelodyTracker::updateMelodyRange(const std::vector<float>& magnitudes, float binWidth) {
    // Calculate new centroid with peak emphasis
    float newCentroid = calculateMelodyCentroid(magnitudes, binWidth);
    float newEnergy = calculateMelodyEnergy(magnitudes, binWidth);

    // Find melody peaks for better note detection
    auto peaks = findMelodyPeaks(magnitudes, binWidth);
    float dominantFreq = peaks.first;

    if (newEnergy > energyThreshold) {
        // Use the dominant peak if it's strong enough
        if (dominantFreq > 0) {
            newCentroid = dominantFreq;
        }

        // Fast adaptation for note changes
        if (currentRange.confidence > minConfidence) {
            currentRange.centroid = currentRange.centroid * (1.0f - adaptationRate) + newCentroid * adaptationRate;
        }
        else {
            // Low confidence, accept new centroid quickly
            currentRange.centroid = currentRange.centroid * 0.3f + newCentroid * 0.7f;
        }

        // Update range boundaries around the centroid
        // Narrow range for better note resolution
        float rangeWidth = 400.0f; // Fixed 400Hz window around dominant frequency
        currentRange.start = std::max(melodyMin, currentRange.centroid - rangeWidth * 0.5f);
        currentRange.end = std::min(melodyMax, currentRange.centroid + rangeWidth * 0.5f);

        // Update energy
        currentRange.energy = currentRange.energy * (1.0f - adaptationRate) + newEnergy * adaptationRate;

        // Increase confidence quickly for strong signals
        currentRange.confidence = std::min(1.0f, currentRange.confidence + adaptationRate * 3.0f);
    }
    else {
        // Low energy, decay confidence
        currentRange.confidence *= confidenceDecay;
        currentRange.energy *= 0.9f;
    }
}

std::vector<float> MelodyTracker::getMelodyBands(const std::vector<float>& magnitudes, float sampleRate, int numBands) {
    std::vector<float> bands(numBands, 0.0f);

    if (magnitudes.empty()) {
        return bands;
    }

    float binWidth = sampleRate / (magnitudes.size() * 2);

    // Determine active melody range
    float rangeStart, rangeEnd;

    if (currentRange.confidence > minConfidence) {
        rangeStart = currentRange.start;
        rangeEnd = currentRange.end;
    }
    else {
        // Fallback to common melody range
        rangeStart = 200.0f;
        rangeEnd = 2000.0f;
    }

    // Ensure minimum useful range
    float minRangeWidth = 300.0f;
    if (rangeEnd - rangeStart < minRangeWidth) {
        float center = (rangeStart + rangeEnd) * 0.5f;
        rangeStart = center - minRangeWidth * 0.5f;
        rangeEnd = center + minRangeWidth * 0.5f;
    }

    // Create high-resolution frequency bands
    float rangeWidth = rangeEnd - rangeStart;
    float bandWidth = rangeWidth / numBands;

    for (int i = 0; i < numBands; i++) {
        float freqStart = rangeStart + i * bandWidth;
        float freqEnd = rangeStart + (i + 1) * bandWidth;

        int startBin = std::max(1, (int)(freqStart / binWidth));
        int endBin = std::min((int)magnitudes.size() - 1, (int)(freqEnd / binWidth));

        if (startBin <= endBin) {
            float energy = 0.0f;
            float peakEnergy = 0.0f;
            int count = 0;

            for (int bin = startBin; bin <= endBin; bin++) {
                float magnitude = magnitudes[bin];

                // Linear aggregation for melody
                energy += magnitude;

                // Detect peaks within this band for note clarity
                if (bin > 0 && bin < magnitudes.size() - 1) {
                    if (magnitude > magnitudes[bin - 1] && magnitude > magnitudes[bin + 1]) {
                        peakEnergy = std::max(peakEnergy, magnitude);
                    }
                }

                count++;
            }

            if (count > 0) {
                float baseEnergy = energy / count;

                // Emphasize bands with peaks (likely notes)
                float finalEnergy = baseEnergy;
                if (peakEnergy > energyThreshold) {
                    finalEnergy = baseEnergy + peakEnergy * 0.5f; // 50% peak boost
                }

                bands[i] = finalEnergy;
            }
        }
    }

    // Light smoothing for stability (less than before for note tracking)
    if (previousBands.size() == numBands) {
        float smoothing = 0.08f; // Very light smoothing
        for (int i = 0; i < numBands; i++) {
            bands[i] = previousBands[i] * (1.0f - smoothing) + bands[i] * smoothing;
        }
    }

    previousBands = bands;
    return bands;
}

std::pair<float, float> MelodyTracker::getCurrentMelodyRange() const {
    if (currentRange.confidence > minConfidence) {
        return { currentRange.start, currentRange.end };
    }
    return { 200.0f, 2000.0f }; // Default range
}