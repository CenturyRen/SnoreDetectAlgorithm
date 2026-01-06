#define _USE_MATH_DEFINES
#include<stdio.h>
#include<vector>
#include<complex>
#include <cmath>
#include <iostream>
#include <fstream>

#include<pocketfft_hdronly.h>
#include "dr_wav.h"


using namespace std;
// 结构体：用来返回音频数据和元信息
struct AudioData {
    std::vector<float> soungAudio; // 音频数据 (已归一化到 -1.0 到 1.0)
    unsigned int sampleRate;
    unsigned int channels;
};

AudioData readWav(const char* filename);
std::vector<float> stereoToMono(const std::vector<float>& stereoData);
vector<float> preEmphasis(vector<float> originSianal);
vector<float> buildMELFilter(int len, int startIndex, int centerIndex, int stopIndex);
std::vector<std::vector<float>> createDCTMatrix(int nFilters, int nCepstral);
std::vector<float> applyDCT(const std::vector<float>& logEnergies, const std::vector<std::vector<float>>& dctMat);
vector<vector<float>> buildMELFilters(int sampleRate, int n_filter);


vector<vector<float>>buildMFCC(AudioData &originSignal,float stepTime,float stepDuration,int n_melFilter,int n_FFT,int numCepstral);
vector<vector<float>>buildMFCC(const vector<float>& sound, int sampleRate = 44100,float stepTime = 0.030 ,float stepDuration = 0.015, int n_melFilter = 32, int n_FFT = 2048, int numCepstral = 13);

// vector<vector<float>>buildMFCC(const std::vector<float>& sound);