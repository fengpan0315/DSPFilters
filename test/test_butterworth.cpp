#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "DspFilters/Dsp.h"

enum FilterType {
    LOWPASS,
    HIGHPASS,
    BANDPASS,
    BANDSTOP
};

struct FilterConfig {
    FilterType type;
    int order;
    double sampleRate;
    double cutoffFrequency;
    double centerFrequency;
    double widthFrequency;
};

void writeBiquadToFile(std::ofstream& file, int stageIndex,
                       double a0, double a1, double a2,
                       double b0, double b1, double b2) {
    file << "Stage " << stageIndex << ": "
         << a0 << " "
         << a1 << " "
         << a2 << " "
         << b0 << " "
         << b1 << " "
         << b2 << "\n";
}

void outputButterworthFilter(const FilterConfig& config, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        printf("Error: Cannot open file %s\n", filename.c_str());
        return;
    }

    printf("\n========================================\n");
    printf("Butterworth Filter Configuration\n");
    printf("========================================\n");

    switch(config.type) {
        case LOWPASS: {
            Dsp::SimpleFilter <Dsp::Butterworth::LowPass <12>, 1> f;
            f.setup(config.order, config.sampleRate, config.cutoffFrequency);

            int numStages = f.getNumStages();
            printf("Number of stages: %d\n\n", numStages);

            file << "Butterworth Lowpass Filter\n";
            file << "Order: " << config.order << "\n";
            file << "Sample Rate: " << config.sampleRate << " Hz\n";
            file << "Cutoff Frequency: " << config.cutoffFrequency << " Hz\n";
            file << "Number of stages: " << numStages << "\n\n";

            for (int i = 0; i < numStages; i++) {
                const auto& stage = f[i];
                writeBiquadToFile(file, i,
                                 stage.getA0(), stage.getA1(), stage.getA2(),
                                 stage.getB0(), stage.getB1(), stage.getB2());

                printf("Stage %d:\n", i);
                printf("  b0=%.10f, b1=%.10f, b2=%.10f\n",
                       stage.getB0(), stage.getB1(), stage.getB2());
                printf("  a0=%.10f, a1=%.10f, a2=%.10f\n\n",
                       stage.getA0(), stage.getA1(), stage.getA2());
            }
            break;
        }

        case HIGHPASS: {
            Dsp::SimpleFilter <Dsp::Butterworth::HighPass <12>, 1> f;
            f.setup(config.order, config.sampleRate, config.cutoffFrequency);

            int numStages = f.getNumStages();
            printf("Number of stages: %d\n\n", numStages);

            file << "Butterworth Highpass Filter\n";
            file << "Order: " << config.order << "\n";
            file << "Sample Rate: " << config.sampleRate << " Hz\n";
            file << "Cutoff Frequency: " << config.cutoffFrequency << " Hz\n";
            file << "Number of stages: " << numStages << "\n\n";

            for (int i = 0; i < numStages; i++) {
                const auto& stage = f[i];
                writeBiquadToFile(file, i,
                                 stage.getA0(), stage.getA1(), stage.getA2(),
                                 stage.getB0(), stage.getB1(), stage.getB2());

                printf("Stage %d:\n", i);
                printf("  b0=%.10f, b1=%.10f, b2=%.10f\n",
                       stage.getB0(), stage.getB1(), stage.getB2());
                printf("  a0=%.10f, a1=%.10f, a2=%.10f\n\n",
                       stage.getA0(), stage.getA1(), stage.getA2());
            }
            break;
        }

        case BANDPASS: {
            Dsp::SimpleFilter <Dsp::Butterworth::BandPass <12>, 1> f;
            f.setup(config.order, config.sampleRate, config.centerFrequency, config.widthFrequency);

            int numStages = f.getNumStages();
            printf("Number of stages: %d\n\n", numStages);

            file << "Butterworth Bandpass Filter\n";
            file << "Order: " << config.order << "\n";
            file << "Sample Rate: " << config.sampleRate << " Hz\n";
            file << "Center Frequency: " << config.centerFrequency << " Hz\n";
            file << "Bandwidth: " << config.widthFrequency << " Hz\n";
            file << "Number of stages: " << numStages << "\n\n";

            for (int i = 0; i < numStages; i++) {
                const auto& stage = f[i];
                writeBiquadToFile(file, i,
                                 stage.getA0(), stage.getA1(), stage.getA2(),
                                 stage.getB0(), stage.getB1(), stage.getB2());

                printf("Stage %d:\n", i);
                printf("  b0=%.10f, b1=%.10f, b2=%.10f\n",
                       stage.getB0(), stage.getB1(), stage.getB2());
                printf("  a0=%.10f, a1=%.10f, a2=%.10f\n\n",
                       stage.getA0(), stage.getA1(), stage.getA2());
            }
            break;
        }

        case BANDSTOP: {
            Dsp::SimpleFilter <Dsp::Butterworth::BandStop <12>, 1> f;
            f.setup(config.order, config.sampleRate, config.centerFrequency, config.widthFrequency);

            int numStages = f.getNumStages();
            printf("Number of stages: %d\n\n", numStages);

            file << "Butterworth Bandstop Filter\n";
            file << "Order: " << config.order << "\n";
            file << "Sample Rate: " << config.sampleRate << " Hz\n";
            file << "Center Frequency: " << config.centerFrequency << " Hz\n";
            file << "Bandwidth: " << config.widthFrequency << " Hz\n";
            file << "Number of stages: " << numStages << "\n\n";

            for (int i = 0; i < numStages; i++) {
                const auto& stage = f[i];
                writeBiquadToFile(file, i,
                                 stage.getA0(), stage.getA1(), stage.getA2(),
                                 stage.getB0(), stage.getB1(), stage.getB2());

                printf("Stage %d:\n", i);
                printf("  b0=%.10f, b1=%.10f, b2=%.10f\n",
                       stage.getB0(), stage.getB1(), stage.getB2());
                printf("  a0=%.10f, a1=%.10f, a2=%.10f\n\n",
                       stage.getA0(), stage.getA1(), stage.getA2());
            }
            break;
        }
    }

    file.close();
    printf("Filter parameters saved to: %s\n", filename.c_str());
    printf("========================================\n");
}

int main()
{
    printf("\n=== Example 1: Butterworth Lowpass Filter ===\n");
    FilterConfig lowpassConfig;
    lowpassConfig.type = LOWPASS;
    lowpassConfig.order = 4;
    lowpassConfig.sampleRate = 44100;
    lowpassConfig.cutoffFrequency = 1000;
    outputButterworthFilter(lowpassConfig, "butterworth_lowpass.txt");

    printf("\n=== Example 2: Butterworth Highpass Filter ===\n");
    FilterConfig highpassConfig;
    highpassConfig.type = HIGHPASS;
    highpassConfig.order = 4;
    highpassConfig.sampleRate = 44100;
    highpassConfig.cutoffFrequency = 500;
    outputButterworthFilter(highpassConfig, "butterworth_highpass.txt");

    printf("\n=== Example 3: Butterworth Bandpass Filter ===\n");
    FilterConfig bandpassConfig;
    bandpassConfig.type = BANDPASS;
    bandpassConfig.order = 4;
    bandpassConfig.sampleRate = 44100;
    bandpassConfig.centerFrequency = 2000;
    bandpassConfig.widthFrequency = 500;
    outputButterworthFilter(bandpassConfig, "butterworth_bandpass.txt");

    printf("\n=== Example 4: Butterworth Bandstop Filter ===\n");
    FilterConfig bandstopConfig;
    bandstopConfig.type = BANDSTOP;
    bandstopConfig.order = 4;
    bandstopConfig.sampleRate = 44100;
    bandstopConfig.centerFrequency = 1000;
    bandstopConfig.widthFrequency = 500;
    outputButterworthFilter(bandstopConfig, "butterworth_bandstop.txt");

    return 0;
}
