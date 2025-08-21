#pragma once
#include "ofMain.h"
#include "Wave.h"
#include "ValueToggler.h"
#include "AudioAnalyzer.h"
#include "AudioVisualizer.h"
#include <iomanip>

class ofApp : public ofBaseApp {
public:
    void setup();
    void update();
    void draw();
    void keyPressed(int key);
    void drawControls();
    void drawVisualizerControls();
    void keyReleased(int key);
    void mouseMoved(int x, int y);
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void mouseEntered(int x, int y);
    void mouseExited(int x, int y);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);

    // Audio input callback (using ofSoundBuffer)
    void audioIn(ofSoundBuffer& input);

    // Audio data storage
    vector<float> left;         // Left channel audio data
    vector<float> right;        // Right channel audio data
    vector<float> volHistory;   // Volume history for visualization

    float leftVol, rightVol;
    float smoothedLeftVol, smoothedRightVol;
    float scaledVol;            // Scaled volume (0-1 range)

    // Audio buffer size
    int bufferSize;             // Audio buffer size


    // Sound stream object
    ofSoundStream soundStream;


    // Wave visualization
    Wave soundWave;
    Wave leftWave;
    Wave rightWave;
    vector<Wave> waves;
    float frequency;
    float amplitude;
    float phase;
    int sampleRate;
    int numSamples;

    float fastSmoothingFactor = 0.1f;   // For quick response
    float slowSmoothingFactor = 0.0002f;  // For stability
    float smoothingThreshold = 0.05f;   // Threshold for switching between modes

    int audioBufferSize;
    float smoothedAmplitude;
    float smoothedFrequency;

    // Mode control
    bool useAudioInput;
    float baseFrequency;
    float baseAmplitude;

    std::unique_ptr<AudioAnalyzer> audioAnalyzer;
    AudioFeatures currentFeatures;

    std::unique_ptr<AudioVisualizer> visualizer;
    VisualizationParams visualParams;
    bool showVisualizerControls = false;
};