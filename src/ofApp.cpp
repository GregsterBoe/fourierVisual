#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup() {
    ofBackground(0);
    ofSetVerticalSync(true);

    center.x = ofGetWindowWidth() / 2;
    center.y = ofGetWindowHeight() / 2;
    topLeft.x = ofGetWindowWidth() / 4;
    topLeft.y = ofGetWindowHeight() / 4;
    topRight.x = ofGetWindowWidth() / 1.5;
    topRight.y = ofGetWindowHeight() / 4;
    bottomLeft.x = ofGetWindowWidth() / 4;
    bottomLeft.y = ofGetWindowHeight() / 1.5;
    bottomRight.x = ofGetWindowWidth() / 1.5;
    bottomRight.y = ofGetWindowHeight() / 1.5;
    midRight.x = 0;
    midRight.y = ofGetWindowHeight() / 2;

    // Initialize wave parameters
    baseFrequency = 440.0; // Base frequency (A4)
    baseAmplitude = 20.0;
    frequency = baseFrequency;
    amplitude = baseAmplitude;
    phase = 0.0;
    sampleRate = 100;
    numSamples = 1024;


    // Setup visualizer
    visualizer = std::make_unique<AudioVisualizer>();
    visualizer->setup(sampleRate, bufferSize);

    // Configure visualization parameters
    visualParams.position = glm::vec2(50, 50);
    visualParams.size = glm::vec2(700, 400);
    visualParams.primaryColor = ofColor::cyan;
    visualParams.secondaryColor = ofColor::magenta;
    visualParams.scale = 1.0f;
    visualParams.smoothing = 0.1f;
    visualParams.sensitivity = 1.0f;
    visualParams.showLabels = true;

    visualizer->setParams(visualParams);

    ofSoundStreamSettings settings;
    // --- Choose your preferred Windows API ---
    // For DirectSound:
    settings.setApi(ofSoundDevice::Api::MS_DS);
    // OR for WASAPI (generally recommended for modern Windows):
    // settings.setApi(ofSoundDevice::Api::MS_WASAPI);


    auto devices = soundStream.getDeviceList(); // Get devices for the chosen API

    int inputDeviceID = -1;
    cout << "Available Audio Devices for chosen API:" << endl;
    for (int i = 0; i < devices.size(); i++) {
        cout << "device: " << devices[i].deviceID << " name: " << devices[i].name
            << " input channels: " << devices[i].inputChannels
            << " output channels: " << devices[i].outputChannels << endl;

        // Look for "VoiceMeeter Output" or "VoiceMeeter Aux Output"
        // These are the *recording* devices you want to listen to.
        if (devices[i].name.find("VoiceMeeter Aux Output (VB-Audio VoiceMeeter AUX VAIO)") != string::npos && devices[i].inputChannels > 0) {
            inputDeviceID = devices[i].deviceID;
            cout << "Found VoiceMeeter Aux Output (WASAPI/DirectSound) at deviceID: " << inputDeviceID << endl;
            break;
        }
        if (devices[i].name.find("VoiceMeeter Output (VB-Audio VoiceMeeter VAIO)") != string::npos && devices[i].inputChannels > 0 && inputDeviceID == -1) {
            inputDeviceID = devices[i].deviceID;
            cout << "Found VoiceMeeter Output (WASAPI/DirectSound) at deviceID: " << inputDeviceID << endl;
            // Don't break yet, in case Aux is found later and preferred
        }
    }

    if (inputDeviceID != -1) {
        settings.setInDevice(devices[inputDeviceID]);
        cout << "Using input device: " << devices[inputDeviceID].name << endl;
    }
    else {
        cout << "VoiceMeeter recording device not found for selected API. Using default.";
    }

    settings.setInListener(this);
    settings.sampleRate = 44100;
#ifdef TARGET_EMSCRIPTEN
    settings.numOutputChannels = 2;
#else
    settings.numOutputChannels = 0;
#endif
    settings.numInputChannels = 2; // Stereo input
    settings.bufferSize = bufferSize;
    soundStream.setup(settings);
    
    // initialize audio analysis
    audioAnalyzer = std::make_unique<AudioAnalyzer>();
    audioAnalyzer->setup(44100, bufferSize);

    // Initialize togglers
    sampleRateToggler = std::make_unique<ValueToggler<int>>(sampleRate, 100, 'S', 's');
    numSamplesToggler = std::make_unique<ValueToggler<int>>(numSamples, 100, 'N', 'n');
    frequencyToggler = std::make_unique<ValueToggler<float>>(baseFrequency, 5.0f, 'F', 'f');
    amplitudeToggler = std::make_unique<ValueToggler<float>>(baseAmplitude, 5.0f, 'A', 'a');
    phaseToggler = std::make_unique<ValueToggler<float>>(phase, 5.0f, 'P', 'p');
}

//--------------------------------------------------------------
void ofApp::audioIn(ofSoundBuffer& input) {
    size_t numFrames = input.getNumFrames();

    if (left.size() != numFrames) {
        left.resize(numFrames);
        right.resize(numFrames);
    }

    // Extract audio data (same as before)
    for (size_t i = 0; i < numFrames && i < left.size(); i++) {
        if (input.getNumChannels() >= 2) {
            left[i] = input.getSample(i, 0);
            right[i] = input.getSample(i, 1);
        }
    }

    // NEW: Analyze audio using current mode
    currentFeatures = audioAnalyzer->analyze(left, right);

    // Update smoothed values (you can now use currentFeatures.leftRMS instead)
    float leftChange = abs(currentFeatures.leftRMS - smoothedLeftVol);
    float rightChange = abs(currentFeatures.rightRMS - smoothedRightVol);

    float leftSmoothing = (leftChange > smoothingThreshold) ? fastSmoothingFactor : slowSmoothingFactor;
    float rightSmoothing = (rightChange > smoothingThreshold) ? fastSmoothingFactor : slowSmoothingFactor;

    smoothedLeftVol = smoothedLeftVol * (1.0f - leftSmoothing) + currentFeatures.leftRMS * leftSmoothing;
    smoothedRightVol = smoothedRightVol * (1.0f - rightSmoothing) + currentFeatures.rightRMS * rightSmoothing;



}

//--------------------------------------------------------------
void ofApp::update() {
    // ... your existing update code ...

    // Update audio analysis
    if (audioAnalyzer && (useAudioInput || !left.empty())) {
        if (useAudioInput) {
            currentFeatures = audioAnalyzer->analyze(left, right);
        }
        else {
            frequency = frequencyToggler->getValue();
            amplitude = amplitudeToggler->getValue();
            leftWave = Wave(frequency, amplitude, phase, sampleRate, numSamples);
            rightWave = Wave(frequency, amplitude, phase, sampleRate, numSamples);
        }

        // Update visualizer with current audio features and raw data
        visualizer->update(currentFeatures, left, right);
    } 

    // Update phase for animation
    phase += 0.1;
}

void ofApp::draw() {
    ofBackground(20);

    // ... your existing drawing code ...

    // Draw the visualizer
    if (visualizer) {
        visualizer->draw();
    }

    // Draw controls
    drawControls();

    if (showVisualizerControls) {
        drawVisualizerControls();
    }
}

void ofApp::drawControls() {
    string info = "Mode: " + string(useAudioInput ? "Audio Input" : "Manual Control");
    info += "\nAnalysis Mode: " + audioAnalyzer->getModeNames()[static_cast<int>(audioAnalyzer->getCurrentMode())];
    info += "\nFrequency: " + ofToString(frequency, 1) + " Hz";

    if (useAudioInput) {
        // Show different information based on analysis mode
        switch (audioAnalyzer->getCurrentMode()) {
        case AudioAnalysisMode::RMS_AMPLITUDE:
            info += "\nLeft RMS: " + ofToString(currentFeatures.leftRMS, 3);
            info += "\nRight RMS: " + ofToString(currentFeatures.rightRMS, 3);
            break;

        case AudioAnalysisMode::FFT_SPECTRUM:
            info += "\nSpectrum Bins: " + ofToString(currentFeatures.fftMagnitudes.size());
            info += "\nOverall Energy: " + ofToString(currentFeatures.overallEnergy, 3);
            break;

        case AudioAnalysisMode::PEAK_DETECTION:
            info += "\nPeak Detected: " + string(currentFeatures.peakDetected ? "YES" : "NO");
            if (currentFeatures.peakDetected) {
                info += "\nPeak Magnitude: " + ofToString(currentFeatures.peakMagnitude, 3);
            }
            break;

        case AudioAnalysisMode::SPECTRAL_CENTROID:
            info += "\nSpectral Centroid: " + ofToString(currentFeatures.spectralCentroid, 1) + " Hz";
            break;

        case AudioAnalysisMode::ONSET_DETECTION:
            info += "\nOnset Detected: " + string(currentFeatures.onsetDetected ? "YES" : "NO");
            if (currentFeatures.onsetDetected) {
                info += "\nOnset Strength: " + ofToString(currentFeatures.onsetStrength, 3);
            }
            break;
        }
    }

    info += "\n\nControls:";
    info += "\n[SPACE] Toggle Audio Input/Manual";
    info += "\n[M/m] Change Analysis Mode";
    info += "\n" + audioAnalyzer->getModeNames()[static_cast<int>(audioAnalyzer->getCurrentMode())];

    ofDrawBitmapString(info, 20, 30);

    // Add visualizer info
    if (visualizer) {
        auto modeNames = visualizer->getModeNames();
        int currentModeIndex = static_cast<int>(visualizer->getCurrentMode());

        ofSetColor(255);
        ofDrawBitmapString("Visualization Mode: " + modeNames[currentModeIndex], 10, ofGetHeight() - 60);
        ofDrawBitmapString("Controls: V=toggle controls, LEFT/RIGHT=change mode, UP/DOWN=sensitivity", 10, ofGetHeight() - 40);
        ofDrawBitmapString("Numbers 1-8: direct mode select, R=reset", 10, ofGetHeight() - 20);
    }
}

void ofApp::drawVisualizerControls() {
    if (!visualizer) return;

    // Draw a semi-transparent overlay
    ofPushStyle();
    ofSetColor(0, 0, 0, 150);
    ofDrawRectangle(ofGetWidth() - 300, 10, 280, 250);

    ofSetColor(255);
    int y = 30;
    ofDrawBitmapString("VISUALIZER CONTROLS", ofGetWidth() - 290, y);
    y += 20;

    // Current mode
    auto modeNames = visualizer->getModeNames();
    int currentModeIndex = static_cast<int>(visualizer->getCurrentMode());
    ofDrawBitmapString("Mode: " + modeNames[currentModeIndex], ofGetWidth() - 290, y);
    y += 15;

    // Parameters
    ofDrawBitmapString("Sensitivity: " + ofToString(visualParams.sensitivity, 2), ofGetWidth() - 290, y);
    y += 15;
    ofDrawBitmapString("Smoothing: " + ofToString(visualParams.smoothing, 2), ofGetWidth() - 290, y);
    y += 15;
    ofDrawBitmapString("Scale: " + ofToString(visualParams.scale, 2), ofGetWidth() - 290, y);
    y += 20;

    // Audio analyzer info
    if (audioAnalyzer) {
        auto analyzerModeNames = audioAnalyzer->getModeNames();
        int analyzerModeIndex = static_cast<int>(audioAnalyzer->getCurrentMode());
        ofDrawBitmapString("Audio Analysis: " + analyzerModeNames[analyzerModeIndex], ofGetWidth() - 290, y);
        y += 15;
    }

    // Instructions
    y += 10;
    ofDrawBitmapString("LEFT/RIGHT: Change mode", ofGetWidth() - 290, y);
    y += 12;
    ofDrawBitmapString("UP/DOWN: Sensitivity", ofGetWidth() - 290, y);
    y += 12;
    ofDrawBitmapString("1-8: Direct mode select", ofGetWidth() - 290, y);
    y += 12;
    ofDrawBitmapString("R: Reset visualizer", ofGetWidth() - 290, y);
    y += 12;
    ofDrawBitmapString("V: Hide this panel", ofGetWidth() - 290, y);

    ofPopStyle();
}

void ofApp::keyPressed(int key) {
    bool handled = false;

    // Toggle between audio input and manual control
    if (key == ' ') {
        useAudioInput = !useAudioInput;
        std::cout << "Switched to " << (useAudioInput ? "Audio Input" : "Manual Control") << " mode" << std::endl;
        handled = true;
    }

    // NEW: Change analysis mode
    if (key == 'M') {
        audioAnalyzer->nextMode();
        handled = true;
    }
    if (key == 'm') {
        audioAnalyzer->previousMode();
        handled = true;
    }

    // Your existing key handling code...
    auto handleToggle = [&](auto& toggler, const std::string& name) {
        if (key == toggler->getIncreaseKey()) {
            toggler->increase();
            std::cout << name << " increased to: " << std::fixed << std::setprecision(2) << toggler->getValue() << std::endl;
            handled = true;
        }
        else if (key == toggler->getDecreaseKey()) {
            toggler->decrease();
            std::cout << name << " decreased to: " << std::fixed << std::setprecision(2) << toggler->getValue() << std::endl;
            handled = true;
        }
    };

    handleToggle(sampleRateToggler, "Sample Rate");
    handleToggle(numSamplesToggler, "Num Samples");
    handleToggle(phaseToggler, "Phase");

    if (!useAudioInput) {
        handleToggle(frequencyToggler, "Frequency");
        handleToggle(amplitudeToggler, "Amplitude");
    }

    if (!handled) {
        std::cout << "Unmapped key: " << static_cast<char>(key) << std::endl;
    }

    // Your existing g/G keys...
    if (key == 'g') {
        soundStream.start();
    }
    if (key == 'G') {
        soundStream.stop();
    }

    switch (key) {
    case 'v': case 'V':
        showVisualizerControls = !showVisualizerControls;
        break;

    case OF_KEY_LEFT:
        if (visualizer) visualizer->previousMode();
        break;

    case OF_KEY_RIGHT:
        if (visualizer) visualizer->nextMode();
        break;

    case OF_KEY_UP:
        visualParams.sensitivity += 0.1f;
        if (visualizer) visualizer->setParams(visualParams);
        break;

    case OF_KEY_DOWN:
        visualParams.sensitivity = std::max(0.1f, visualParams.sensitivity - 0.1f);
        if (visualizer) visualizer->setParams(visualParams);
        break;

    case 'r': case 'R':
        if (visualizer) visualizer->reset();
        break;

        // Cycle through visualization modes with number keys
    case '1': if (visualizer) visualizer->setMode(VisualizationMode::WAVE_RMS); break;
    case '2': if (visualizer) visualizer->setMode(VisualizationMode::SPECTRUM_BARS); break;
    case '3': if (visualizer) visualizer->setMode(VisualizationMode::CIRCULAR_SPECTRUM); break;
    case '4': if (visualizer) visualizer->setMode(VisualizationMode::PEAK_PULSES); break;
    case '5': if (visualizer) visualizer->setMode(VisualizationMode::CENTROID_WAVE); break;
    case '6': if (visualizer) visualizer->setMode(VisualizationMode::ONSET_PARTICLES); break;
    case '7': if (visualizer) visualizer->setMode(VisualizationMode::WAVEFORM_RAW); break;
    case '8': if (visualizer) visualizer->setMode(VisualizationMode::COMBINED_VIEW); break;
    }
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {
    waves.clear();
    soundWave = Wave(frequency, amplitude, phase, sampleRate, numSamples);
    waves.push_back(soundWave);
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {
}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {
}