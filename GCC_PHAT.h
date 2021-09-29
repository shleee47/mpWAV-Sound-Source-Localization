#ifndef _GCC_PHAT_
#define _GCC_PHAT_
/*
 * Ver 1.0 2021-01-25 init
 */
#define _USE_MATH_DEFINES // to use M_PI
#include <math.h>
#include <tgmath.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex>
#include "STFT/cpp/STFT.h"
#include "STFT/cpp/Ooura_FFT.h"

class GCC_PHAT {
private:
    const double eps = 2.2204e-16;
    const double sound_speed = 343.3;
    const double smooth_factor = 0.50;
    const double IntpRatio = 16.0;
    const int corr_opt = 1;
    const int phat_opt = 1;
    const int interest_freq_opt = 0;
    const double pi = 3.1415926535897;
    
    int ch;
    int fs;
    int nfft;
    int nhfft;
    int shift;
    int frame;
    double max_idx;
    double max_delay;
    double max;
    double mic_dist;

    std::complex<double>* R_tmp;
    std::complex<double>* R_tmp1;
    std::complex<double>* R_tmp2;
    double* R_PHAT_tmp3;
    double* R_PHAT;
    double* R_f_tmp;
    double* R_f_input;
    double* R_f;

    Ooura_FFT* fft = nullptr;
    std::complex<double> i1;

public:

    inline GCC_PHAT(int ch, int nfft, int fs, int shift, double mic_dist);
    inline ~GCC_PHAT();
    inline void Process(double**);

    double Estimated_Sample_delay;
    std::complex<double> Estimated_azimuth;
    double* R_f_output;
};

GCC_PHAT::GCC_PHAT(
    int ch_,
    int nfft_,
    int fs_,
    int shift_,
    double  mic_dist_):i1(0,1){

    ch = ch_;
    nfft = nfft_;
    fs = fs_;
    shift = shift_;
    frame = nfft_;
    mic_dist = mic_dist_;
    
    nhfft = nfft / 2 + 1;
    max_delay = ceil(mic_dist * double(fs) / sound_speed);
    
    R_tmp = new std::complex<double>[nhfft];
    R_tmp1 = new std::complex<double>[nhfft];
    R_tmp2 = new std::complex<double>[nfft / 2 * (IntpRatio)+1];
    R_PHAT_tmp3 = new double[nfft * IntpRatio + 2]; // for transformation right before iFFT
    R_PHAT = new double[nfft * IntpRatio];
    R_f_tmp = new double[max_delay * IntpRatio * 2 + 1];
    R_f_input = new double [(max_delay * IntpRatio) * 2 + 1]; // the very before frame
    R_f = new double[(max_delay * IntpRatio) * 2 + 1];
    R_f_output = new double[(max_delay * IntpRatio) * 2 + 1];
    memset(R_f_output, 0, sizeof(double) * ((max_delay * IntpRatio) * 2 + 1));
    fft = new Ooura_FFT(nfft * IntpRatio, ch);
}

void GCC_PHAT::Process(double** X) {
    int i, j, k, l;
    double weight_PHAT;

    /* Cross-Correlation */
    for (j = 0; j < nhfft; j++) {
        R_tmp[j] = (X[0][j + j] * X[1][j + j]) + (X[0][j + j + 1] * X[1][j + j + 1]) + i1*((X[0][j + j + 1] * X[1][j + j]) - (X[0][j + j ] * X[1][j + j + 1]));   // (ac+bd) + (bc -ad)i
    }

    /* Interest Frequency Limit */
    if (interest_freq_opt == 1) {
        R_tmp[0] = 0;
        for (i = 99; i < nhfft; i++) {
            R_tmp[i] = 0;
        }
    }

    /* PHAT */
    if (phat_opt == 1) {
        for (j = 0; j < nhfft; j++) {
            weight_PHAT = sqrt((R_tmp[j].real() * R_tmp[j].real()) + (R_tmp[j].imag() * R_tmp[j].imag()))+eps;
            R_tmp1[j] = R_tmp[j] / weight_PHAT;
        }
    }
    else {
        for (j = 0; j < nhfft; j++) {
            R_tmp1[j] = R_tmp[j];
        }
    }

    /* Interpolation (Zero Padding) */
    for (j = 0; j < nfft/2 ; j++) {
        R_tmp2[j] = R_tmp1[j];
    }
    R_tmp2[nfft/2* int(IntpRatio)] = R_tmp1[nfft/2]; //중앙에 값 하나

    /* Transform double complex to double*/
    for (j = 0; j < (nfft/2 * IntpRatio)+1; j++) {
        R_PHAT_tmp3[j + j] = R_tmp2[j].real();
        R_PHAT_tmp3[j + j + 1] = R_tmp2[j].imag();
    }

    /*Inverse FFT*/
    fft->SingleiFFT(R_PHAT_tmp3);
    
    ///*Half Flip*/  마지막 두개는 버려
    for (j = 0; j < nfft * IntpRatio / 2; j++) {
        R_PHAT[nfft * int(IntpRatio) / 2 + j] = R_PHAT_tmp3[j];
        R_PHAT[j] = R_PHAT_tmp3[nfft * int(IntpRatio) / 2 + j];
    }

    ///* Sample Delay Estimation */
    for (j = 0; j < max_delay * IntpRatio * 2 + 1; j++) {
        R_f_tmp[j] = R_PHAT[(nfft * int(IntpRatio) / 2) - 1 - int(max_delay * IntpRatio) + j];
    }

    /* Smoothing Factor */
    if (corr_opt == 1) {
        for (j = 0; j < (max_delay * IntpRatio) * 2 + 1; j++) {
            R_f[j] = R_f_tmp[j];
        }
    }
    else if (corr_opt == 2) {
        for (j = 0; j < (max_delay * IntpRatio) * 2 + 1; j++) {
            R_f[j] = smooth_factor * R_f_output[j] + (1 - smooth_factor) * R_f_tmp[j];
        }
    }
    else if (corr_opt == 3) {
        for (j = 0; j < (max_delay * IntpRatio) * 2 + 1; j++) {
            R_f[j] = R_f_output[j]+ R_f_tmp[j];
        }
    }

    /* Sample Delay Estimation */
    max = R_f[0];
    max_idx = 0;
    for (j = 1; j < (max_delay * IntpRatio) * 2 + 1 ; j++) {
        if (R_f[j] > max) {
            max_idx = j;
            max = R_f[j];
        }
    }

    /* Public Variable for return */
    Estimated_Sample_delay = (max_idx) / IntpRatio - max_delay;
    //Estimated_azimuth = (asin(Estimated_Sample_delay * sound_speed / double(fs) / mic_dist)) * 180 / pi;
    //printf("(max_delay*IntpRatio)*2 : %f\n", (max_delay*IntpRatio)*2);
    //printf("max_delay : %f\n", max_delay);
    //printf("Estimated_Sample_delay : %f\n", Estimated_Sample_delay);
    Estimated_azimuth = (asin(Estimated_Sample_delay / max_delay)) * 180 / pi;
    for (j = 0; j < (max_delay * IntpRatio) * 2 + 1; j++) {
        R_f_output[j] = R_f[j];
    }
}

/* Delete Function */
GCC_PHAT::~GCC_PHAT() {
    int i, j, k, l;
    delete[] R_f_input;
    delete[] R_f;
    delete[] R_tmp;
    delete[] R_tmp1;
    delete[] R_tmp2;
    delete[] R_PHAT;
    delete[] R_f;
    delete[] R_PHAT_tmp3;
}
#endif
