/*
  ==============================================================================

    This file was auto-generated by the Jucer!

    It contains the basic startup code for a Juce application.

  ==============================================================================
*/

#ifndef __PLUGINPROCESSOR_H_4FD11FCE__
#define __PLUGINPROCESSOR_H_4FD11FCE__

#include "JuceHeader.h"
#include "JucePluginCharacteristics.h"


//==============================================================================
/**
*/
class TheFunctionAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    TheFunctionAudioProcessor();
    ~TheFunctionAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

    void processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages) override;

    //==============================================================================
#if ! JUCE_AUDIOPROCESSOR_NO_GUI
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;
#endif

    double getTailLengthSeconds() const override { return 0.0; }

    //==============================================================================
    const String getName() const override;

    int getNumParameters() override;

    float getParameter (int index) override;
    void setParameter (int index, float newValue) override;

    const String getParameterName (int index) override;
    const String getParameterText (int index) override;

    const String getInputChannelName (int channelIndex) const override;
    const String getOutputChannelName (int channelIndex) const override;
    bool isInputChannelStereoPair (int index) const override;
    bool isOutputChannelStereoPair (int index) const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    bool silenceInProducesSilenceOut() const override { return true; }

    enum Parameters
    {
        gainParam = 0,
        gainLParam,
        gainRParam,
        phaseLParam,
        phaseRParam,
        panLParam,
        panRParam,
        totalNumParams
    };

    float gain;
    float gainL;
    float gainR;

    float panL;
    float panR;

    float phaseL;
    float phaseR;

    int currentPreset;
    int timeSinceChunkCalled;

    AudioSampleBuffer tmpBuffer;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TheFunctionAudioProcessor);
};

#endif  // __PLUGINPROCESSOR_H_4FD11FCE__