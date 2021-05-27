// clang-format off

/*
	==============================================================================
	This file is part of Tal-Reverb by Patrick Kunz.

	Copyright(c) 2005-2009 Patrick Kunz, TAL
	Togu Audio Line, Inc.
	http://kunz.corrupt.ch

	This file may be licensed under the terms of of the
	GNU General Public License Version 2 (the ``GPL'').

	Software distributed under the License is distributed
	on an ``AS IS'' basis, WITHOUT WARRANTY OF ANY KIND, either
	express or implied. See the GPL for the specific language
	governing rights and limitations.

	You should have received a copy of the GPL along with this
	program. If not, go to http://www.gnu.org/licenses/gpl.html
	or write to the Free Software Foundation, Inc.,
	51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
	==============================================================================
 */

// clang-format off
#if !defined(__TalReverb_h)
#define __TalReverb_h

#include <optional>
#include <array>
#include <vector>

#include "AllPassFilter.h"
#include "CombFilter.h"
#include "NoiseGenerator.h"
#include "Filter.h"
#include "math.h"
#include "AudioUtils.h"
#include "TalEq.h"

namespace artv_dsp_pull { namespace tal_reverb2 {

class TalReverb
{
private:
	static constexpr int DELAY_LINES_COMB = 4;
	static constexpr int DELAY_LINES_ALLPASS = 5;

	static constexpr int MAX_PRE_DELAY_MS = 1000;
	std::vector<float> reflectionGains;
	std::vector<float> reflectionDelays;

	std::optional<CombFilter> combFiltersPreDelayL;
	std::optional<CombFilter> combFiltersPreDelayR;

	std::array<std::optional<CombFilter>, DELAY_LINES_COMB> combFiltersL;
	std::array<std::optional<CombFilter>,DELAY_LINES_COMB> combFiltersR;
	std::array<std::optional<NoiseGenerator>, DELAY_LINES_COMB> noiseGeneratorAllPassL;
	std::array<std::optional<NoiseGenerator>, DELAY_LINES_COMB> noiseGeneratorAllPassR;
	std::array<std::optional<NoiseGenerator>, DELAY_LINES_COMB> noiseGeneratorDelayL;
	std::array<std::optional<NoiseGenerator>, DELAY_LINES_COMB> noiseGeneratorDelayR;
	std::array<std::optional<NoiseGenerator>, DELAY_LINES_COMB> diffusionL;
	std::array<std::optional<NoiseGenerator>, DELAY_LINES_COMB> diffusionR;

	std::array<std::optional<AllPassFilter>, DELAY_LINES_ALLPASS> allPassFiltersL;
	std::array<std::optional<AllPassFilter>, DELAY_LINES_ALLPASS> allPassFiltersR;

	std::optional<AllPassFilter> preAllPassFilterL;
	std::optional<AllPassFilter> preAllPassFilterR;

	std::optional<AllPassFilter> postAllPassFilterL;
	std::optional<AllPassFilter> postAllPassFilterR;

	std::optional<TalEq> talEqL;
	std::optional<TalEq> talEqR;

	float decayTime;
	float preDelayTime;
	bool stereoMode;
	float modulationIntensity;

	float outL;
	float outR;

	AudioUtils audioUtils;

public:
	TalReverb(int sampleRate)
	{
		createDelaysAndCoefficients(DELAY_LINES_COMB + DELAY_LINES_ALLPASS, 82.0f);

		combFiltersPreDelayL.emplace((float)MAX_PRE_DELAY_MS, 0.0f, sampleRate);
		combFiltersPreDelayR.emplace((float)MAX_PRE_DELAY_MS, 0.0f, sampleRate);

		float stereoSpreadValue = 0.008f;
		float stereoSpreadSign = 1.0f;
		for (int i = 0; i < DELAY_LINES_COMB; i++)
		{
			float stereoSpreadFactor = 1.0f + stereoSpreadValue;
			if (stereoSpreadSign > 0.0f)
			{
				combFiltersL[i].emplace(reflectionDelays[i] * stereoSpreadFactor, reflectionGains[i], sampleRate);
				combFiltersR[i].emplace(reflectionDelays[i], reflectionGains[i], sampleRate);
			}
			else
			{
				combFiltersL[i].emplace(reflectionDelays[i], reflectionGains[i], sampleRate);
				combFiltersR[i].emplace(reflectionDelays[i] * stereoSpreadFactor, reflectionGains[i], sampleRate);
			}
			stereoSpreadSign *= -1.0f;
			noiseGeneratorAllPassL[i].emplace(sampleRate);
			noiseGeneratorAllPassR[i].emplace(sampleRate);
			noiseGeneratorDelayL[i].emplace(sampleRate);
			noiseGeneratorDelayR[i].emplace(sampleRate);
			diffusionL[i].emplace(sampleRate);
			diffusionR[i].emplace(sampleRate);
		}

		preAllPassFilterL.emplace(20.0f,  0.68f, sampleRate);
		preAllPassFilterR.emplace(20.0f,  0.68f, sampleRate);

		postAllPassFilterL.emplace(200.0f,  0.68f, sampleRate);
		postAllPassFilterR.emplace(200.0f,  0.68f, sampleRate);

		for (int i = 0; i < DELAY_LINES_ALLPASS; i++)
		{
			allPassFiltersL[i].emplace(reflectionDelays[i + DELAY_LINES_COMB - 1] * 0.105f,  0.68f, sampleRate);
			allPassFiltersR[i].emplace(reflectionDelays[i + DELAY_LINES_COMB - 1] * 0.1f,  0.68f, sampleRate);
        }

		talEqL.emplace(sampleRate);
		talEqR.emplace(sampleRate);

		decayTime = 0.5f;
		preDelayTime = 0.0f;
		modulationIntensity = 0.12f;
		stereoMode = false;

		outL = 0.0f;
		outR = 0.0f;
	}

	~TalReverb()
	{}

	void setDecayTime(float decayTime)
	{
		this->decayTime = audioUtils.getLogScaledValueInverted(decayTime) * 0.99f;
	}

#ifdef __MOD_DEVICES__
	void setPreDelay(float preDelayTime)
	{
		this->preDelayTime = preDelayTime * 0.001;
	}
#else
	void setPreDelay(float preDelayTime)
	{
		this->preDelayTime = audioUtils.getLogScaledValue(preDelayTime);
	}
#endif
	void setStereoMode(bool stereoMode)
	{
		this->stereoMode = stereoMode;
	}

	void setLowShelfGain(float lowShelfGain)
	{
		talEqL->setLowShelfGain(lowShelfGain);
		talEqR->setLowShelfGain(lowShelfGain);
	}

	void setHighShelfGain(float highShelfGain)
	{
		talEqL->setHighShelfGain(highShelfGain);
		talEqR->setHighShelfGain(highShelfGain);
	}

	void setLowShelfFrequency(float lowShelfFrequency)
	{
		talEqL->setLowShelfFrequency(lowShelfFrequency);
		talEqR->setLowShelfFrequency(lowShelfFrequency);
	}

	void setHighShelfFrequency(float highShelfFrequency)
	{
		talEqL->setHighShelfFrequency(highShelfFrequency);
		talEqR->setHighShelfFrequency(highShelfFrequency);
	}

	void setPeakFrequency(float peakFrequency)
	{
		talEqL->setPeakFrequency(peakFrequency);
		talEqR->setPeakFrequency(peakFrequency);
	}

	void setPeakGain(float peakGain)
	{
		talEqL->setPeakGain(peakGain);
		talEqR->setPeakGain(peakGain);
	}

	// All input values [0..1]
	inline void process(float* sampleL, float* sampleR)
	{
		float revR;
		float revL;

		if (!stereoMode)
		{
			revL = (*sampleL + *sampleR) * 0.25f;
			revL = combFiltersPreDelayL->process(revL, 0.0f, 0.0f, preDelayTime);
			talEqL->process(&revL);
			revR = revL;
		}
		else
		{
			revL = combFiltersPreDelayL->process(*sampleL * 0.5f, 0.0f, 0.0f, preDelayTime);
			revR = combFiltersPreDelayL->process(*sampleR * 0.5f, 0.0f, 0.0f, preDelayTime);
			talEqL->process(&revL);
			talEqR->process(&revR);
		}

		// ----------------- Comb Filter --------------------
		outL = 0.0f;
		outR = 0.0f;

		float scaledRoomSize = decayTime * 0.998f;
		float sign = 1.0f;
		for (int i = 0; i < DELAY_LINES_COMB; i++)
		{
            outL += sign * combFiltersL[i]->processInterpolated(revL, diffusionL[i]->tickFilteredNoiseFast() * 0.2f, scaledRoomSize, scaledRoomSize + 0.012f * noiseGeneratorDelayL[i]->tickFilteredNoise());
			outR += sign * combFiltersR[i]->processInterpolated(revR, diffusionR[i]->tickFilteredNoiseFast() * 0.2f, scaledRoomSize, scaledRoomSize + 0.012f * noiseGeneratorDelayR[i]->tickFilteredNoise());
			sign *= -1.0f;
		}

		// ----------------- Pre AllPass --------------------
		outL += 0.5f * preAllPassFilterL->processInterpolated(revR, 0.8f + 0.2f * noiseGeneratorAllPassL[0]->tickFilteredNoise(), 0.69f, true);
		outR += 0.5f * preAllPassFilterR->processInterpolated(revL, 0.8f + 0.2f * noiseGeneratorAllPassR[0]->tickFilteredNoise(), 0.69f, true);

		//// ----------------- Post AllPass --------------------
		outL += 0.45f * postAllPassFilterL->processInterpolated(revL, 0.8f + 0.2f * noiseGeneratorAllPassL[1]->tickFilteredNoise(), 0.69f, false);
		outR += 0.45f * postAllPassFilterR->processInterpolated(revR, 0.8f + 0.2f * noiseGeneratorAllPassR[1]->tickFilteredNoise(), 0.69f, false);

		// ----------------- AllPass Filter ------------------
		for (int i = 0; i < DELAY_LINES_ALLPASS; i++)
		{
			outL = allPassFiltersL[i]->process(outL);
			outR = allPassFiltersR[i]->process(outR);
		}

		// ----------------- Write to output / Stereo --------
		*sampleL = outL;
		*sampleR = outR;
	}

	void createDelaysAndCoefficients(int numlines, float delayLength)
	{
		reflectionDelays.clear();
		reflectionDelays.resize(numlines);
		reflectionGains.clear();
		reflectionGains.resize(numlines);

		float volumeScale = (float)(-3.0 * delayLength / log10f(0.5f));
		for (int n = numlines - 1; n >= 0; n--)
		{
			reflectionDelays[numlines -1 - n] = delayLength / powf(2.0f, (float)n / numlines);
			reflectionGains[numlines -1 - n] = powf(10.0f, - (3.0f * reflectionDelays[numlines -1 - n]) / volumeScale);
		}
	}
};

}}

#endif
