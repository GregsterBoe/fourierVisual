#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup() {
    ofBackground(0);
    ofSetVerticalSync(true);

    // Initialize wave parameters
    baseFrequency = 440.0; // Base frequency (A4)
    baseAmplitude = 20.0;
    frequency = baseFrequency;
    amplitude = baseAmplitude;
    phase = 0.0;
    sampleRate = 44100;
    numSamples = 1024;
    bufferSize = 2048;


    // Setup visualizer
    visualizer = std::make_unique<AudioVisualizer>();
    visualizer->setup(sampleRate, bufferSize);

    // Configure visualization parameters
    visualParams.position = glm::vec2(200, 200);
    visualParams.size = glm::vec2(ofGetWindowWidth(), ofGetWindowHeight() / 2);
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


    settings.setInListener(this);
    settings.sampleRate = sampleRate;
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
    audioAnalyzer->setup(sampleRate, bufferSize);
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

}

//--------------------------------------------------------------
void ofApp::update() {
    // ... your existing update code ...

    // Update audio analysis
    if (!left.empty()) {
        currentFeatures = audioAnalyzer->analyze(left, right);
        // Update visualizer with current audio features and raw data
        visualizer->update(currentFeatures);
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
    string info = "";
    info += "Avg Freq: " + std::to_string(audioAnalyzer->getAverageFrequency());
    info += ", Max Freq: " + std::to_string(audioAnalyzer->getMaxFrequency());
    ofSetColor(255);

    ofDrawBitmapString("Analysis info: " + info, 10, ofGetHeight() - 80);


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

    // Instructions
    y += 10;
    ofDrawBitmapString("LEFT/RIGHT: Change mode", ofGetWidth() - 290, y);
    y += 12;
    ofDrawBitmapString("UP/DOWN: Sensitivity", ofGetWidth() - 290, y);
    y += 12;
    ofDrawBitmapString("1-8: Direct mode select", ofGetWidth() - 290, y);
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
    case '7': if (visualizer) visualizer->setMode(VisualizationMode::COMBINED_VIEW); break; 
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