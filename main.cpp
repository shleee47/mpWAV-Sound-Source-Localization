#define _CRT_SECURE_NO_WARNINGS
#include "RT_Input.h"
#include "WAV/WAV.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "STFT/cpp/STFT.h"
#include "SRP_DSB.h"
#include "GCC_PHAT.h"
#include <math.h>

#define PI 3.1415926535897
#define DEVICE 4
#define CHANNEL 2
#define SAMPLE_RATE 16000
#define FRAME_SIZE 512
#define SHIFT_SIZE 128
#define MIC_DISTANCE 0.05

double getRadian(int _num) { return _num * (PI / 180); }

int main() {

    /* Available Device Information */
    RtAudio audio;
    unsigned int devices = audio.getDeviceCount();
    RtAudio::DeviceInfo info;
    for (unsigned int i = 0; i < devices; i++) {
        info = audio.getDeviceInfo(i);
        std::cout << "device = " << i;
        std::cout << ": device name = " << info.name << "\n";
    }

    /* Declaration */
    const int device = DEVICE;
    const int ch = CHANNEL;
    const int fs = SAMPLE_RATE;
    const int frame = FRAME_SIZE;
    const int shift = SHIFT_SIZE;
    const double mic_dist = MIC_DISTANCE;

    /* Initialization */
    RT_Input input(device, 16, fs, shift, frame);
    GCC_PHAT* gcc = nullptr;
    STFT* stft = nullptr;
    WAV* wav_1 = nullptr;
    WAV* wav_2 = nullptr;

    double sample_delay;
    double** raw;
    double** raw_tmp;
    double** data;

    stft = new STFT(ch, frame, shift);
    gcc = new GCC_PHAT(ch, frame, fs, shift, mic_dist);

    raw_tmp = new double* [16];
    for (int i = 0; i < 16; i++)
        raw_tmp[i] = new double[shift];

    raw = new double* [ch];
    for (int i = 0; i < ch; i++)
        raw[i] = new double[shift];

    data = new double* [ch];
    for (int i = 0; i < ch; i++)
        data[i] = new double[frame + 2];

    /* process */
    int unit = shift * 1;
    size_t read_1, read_2;
    bool prev_sil = true;
    bool sil = true;
    int t = 0;

    input.Start();

    while (input.IsRunning()) {
        if (input.data.stock.load() >= shift) {
            input.Convert2Array(raw_tmp);
            for (int i = 0; i < shift; i++) {
                raw[0][i] = raw_tmp[4][i];
                raw[1][i] = raw_tmp[8][i];
            }

            stft->stft(raw, data);
            /* Estimation Result */
            if (!stft->sil) {
                gcc->Process(data);
                printf("t : %d |     Estimated Azimuth %f\n", t, (gcc->Estimated_azimuth.real() + 90));
            }
            if (!prev_sil && stft->sil)
                printf("silence\n");
            prev_sil = stft->sil;
            t++;
        }
        else {
            SLEEP(10);
        }
    }
    printf("Streaming Finished!!!\n");
    delete gcc;
    delete stft;

    for (int i = 0; i < ch; i++) {
        delete[] raw[i];
        delete[] data[i];
    }
    delete[] raw;
    delete[] data;

    input.Stop();
    return 0;
}