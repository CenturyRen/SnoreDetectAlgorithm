#define _USE_MATH_DEFINES
#define DR_WAV_IMPLEMENTATION 

#include<stdio.h>
#include<pocketfft_hdronly.h>
#include<vector>
#include<complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include "dr_wav.h"
using namespace std;

// 结构体：用来返回音频数据和元信息
struct AudioData {
    vector<float> soungAudio; // 音频数据 (已归一化到 -1.0 到 1.0)
    unsigned int sampleRate;
    unsigned int channels;
};


// 函数：读取 WAV 文件
AudioData readWav(const char* filename) {
    unsigned int channels;
    unsigned int sampleRate;
    drwav_uint64 totalPCMFrameCount;

    // 关键函数：直接读取并转换为 float
    // pSampleData 是库分配的内存，用完需要 free
    float* pSampleData = drwav_open_file_and_read_pcm_frames_f32(
        filename, 
        &channels, 
        &sampleRate, 
        &totalPCMFrameCount, 
        NULL
    );

    if (pSampleData == NULL) {
        cerr << "错误: 无法打开或读取文件: " << filename << endl;
        return {}; // 返回空
    }

    // 计算总样本数 = 帧数 * 通道数
    size_t totalSamples = (size_t)(totalPCMFrameCount * channels);

    // 拷贝数据到 vector 中 (方便 C++ 管理)
    AudioData data;
    data.soungAudio.assign(pSampleData, pSampleData + totalSamples);
    data.sampleRate = sampleRate;
    data.channels = channels;

    // 释放库分配的内存
    drwav_free(pSampleData, NULL);

    return data;
}

// 将立体声 (L, R, L, R...) 转为 单声道 (M, M...)
std::vector<float> stereoToMono(const std::vector<float>& stereoData) {
    // 1. 检查数据是否为空
    if (stereoData.empty()) return {};

    // 2. 计算单声道的长度 (总长度除以通道数 2)
    size_t numFrames = stereoData.size() / 2;
    
    // 3. 预分配内存 (避免 push_back 带来的多次内存重分配，提高性能)
    std::vector<float> monoData(numFrames);

    // 4. 遍历并计算平均值
    for (size_t i = 0; i < numFrames; ++i) {
        // i * 2     是左声道 (Left)
        // i * 2 + 1 是右声道 (Right)
        float left  = stereoData[i * 2];
        float right = stereoData[i * 2 + 1];

        // 乘 0.5f 比除以 2.0f 稍微快一点点 (在某些 CPU 上)
        monoData[i] = (left + right) * 0.5f;
    }

    return monoData;
}

vector<float> preEmphasis(vector<float> originSianal)
{
    if (originSianal.empty()) return {};
    vector<float> afterPrcSiganl;
    
    afterPrcSiganl.reserve(originSianal.size());

    afterPrcSiganl.push_back(originSianal.at(0));

    for(int i = 1; i < originSianal.size() ; ++i)
    {
        afterPrcSiganl.push_back(originSianal.at(i) - 0.95 * originSianal.at(i - 1));
    }
    return afterPrcSiganl;
}

vector<float> buildMELFilter(int len, int startIndex, int centerIndex, int stopIndex)
{
    if(startIndex < 0 || centerIndex < 0)return {};
    vector<float> melFilter(len);
    std::fill(melFilter.begin(), melFilter.end(), 0);
    // y = h * (x - startIndex) / (centerIndex - startIndex)
    // height = 2 / (startIndex - stopIndex);
    // y1 = 2*(x-startIndex)/(stopIndex-startIndex)*(centerIndex-startIndex);
    melFilter[centerIndex] = 2.0f /  (stopIndex - startIndex);
    for(int i = startIndex; i < centerIndex; ++i)
    {
        melFilter[i] = 2.0f*(i - startIndex)/((stopIndex - startIndex)*(centerIndex - startIndex));
    }
    for(int i = centerIndex + 1; i < stopIndex; ++i)
    {
        melFilter[i] = 2.0f*(i - stopIndex)/((stopIndex - startIndex)*(centerIndex - stopIndex));
    }
    return melFilter;
}
vector<vector<float>> buildMELFilters(int sampleRate, int n_filter)
{
    if(sampleRate < 0 || n_filter < 0) return {};
    vector<vector<float>> melFilters;
    // mel(f) = 2959*log10(1+f/700)
    // 计算mel域范围
    float startFreq = 2595 * log10(1 + 0/700);
    float stopFreq  = 2595 * log10(1 + sampleRate/(2*700));

    // 
    float step = 1.0f * (stopFreq - startFreq) / (n_filter + 1);

    vector<float>point;
    vector<int>lastPoint;
    for(int i = 0; i <= n_filter + 1; ++i)
    {
        point.push_back(i * step);
    }

    // mel转hz
    // f = 700*(10^(melf/2595)-1)
    for(int i = 0; i < point.size(); ++i)
    {
        lastPoint.push_back(((700*(( pow(10,(point[i])/2595)) - 1) * 2048)/44100));
    }

    int startIndex;
    int centerIndex;
    int stopIndex;

    for(int i = 1; i < n_filter+1; ++i)
    {
        startIndex  = lastPoint[i-1];
        centerIndex = lastPoint[i];
        stopIndex   = lastPoint[i+1];

        melFilters.push_back(buildMELFilter(1025, startIndex, centerIndex, stopIndex));

    }
    
    
    return melFilters;
}


void saveToBinary(const std::vector<std::vector<float>>& data, const std::string& filename) {
    if (data.empty()) return;

    std::ofstream outFile(filename, std::ios::binary);
    
    int rows = (int)data.size();
    int cols = (int)data[0].size();

    // 1. 先写入头部信息：行数和列数
    outFile.write((char*)&rows, sizeof(int));
    outFile.write((char*)&cols, sizeof(int));

    // 2. 逐行写入数据
    for (const auto& row : data) {
        // 直接将 vector 内存块写入文件
        outFile.write((char*)row.data(), cols * sizeof(float));
    }

    outFile.close();
    std::cout << "save seccess:" << filename << std::endl;
}



vector<vector<float>>buildMFCC(AudioData originSignal,)
{

}

int main()
{

    // 文件加载
    const char* filepath = "0_0.wav"; 
    printf("loading: %s ...\n", filepath);
    AudioData originSignal = readWav(filepath);
    if (originSignal.soungAudio.empty()) {
        printf("error 1");
        getchar();
        return -1;
    }
    cout << "load success." << endl;


    
    // 双声道转单声道
    vector<float> momoSignal = stereoToMono(originSignal.soungAudio);
    
    // 定义MFCC信息
    unsigned int sampleRate = originSignal.sampleRate;

    float stepTime     = 0.030;// 30ms
    float stepDuration = 0.015;// 15ms重叠
    // int n_MFCC = int(1.0f * momoSignal.size()/(sampleRate*stepTime));// 1秒以30毫秒均分，特征数设为32
    
    int n_melFilter = 32;
    int n_FFT  = 2048;// FFT系数数量
    int win_Length = int(stepTime * sampleRate);// 窗口宽度（30ms的宽度）
    // int stepLength = int(momoSignal.size() / n_MFCC );
    int stepLength = int(stepDuration * sampleRate);
    int n_MFCC = (momoSignal.size() - win_Length) / stepLength + 1;
    // 预加强
    vector<float> preE_Signal = preEmphasis(momoSignal);

    // 分帧
    vector<vector<float>> frameSplitSignal;
    frameSplitSignal.reserve(n_MFCC); // 预留外层空间，提高效率

    // saveToBinary(preE_Signal,"preE_Signal.bin");

    auto startIndex = preE_Signal.begin();
    auto endIndex = startIndex;
    for(int i = 0; i < n_MFCC; ++i)
    {
        endIndex = startIndex + win_Length;
        frameSplitSignal.emplace_back(startIndex, endIndex);
        startIndex += stepLength;
    }

    // 创建汉明窗，同时乘原始数据
    vector<float> hammingWindow;
    for(int i = 0; i < win_Length; ++i)
    {
        float hwValue = 0.54 - 0.46 * cos(2 * M_PI * i / (win_Length - 1));
        hammingWindow.push_back(hwValue);
    }
    for(int i = 0; i < n_MFCC; ++i)
    {
        for(int j = 0; j < win_Length; ++j)
        {
            frameSplitSignal.at(i)[j] *= hammingWindow[j];
        }
    }

    // auto maxIt = std::max_element(frameSplitSignal[10].begin(), frameSplitSignal[10].end());
    // auto minIt = std::min_element(frameSplitSignal[10].begin(), frameSplitSignal[10].end());
    // cout << "max:" << *maxIt <<  " " << *minIt << endl;
    cout << "save origin data." << endl;
    saveToBinary(frameSplitSignal,"originData.bin");


    size_t spectrumSize = n_FFT / 2 + 1;
    vector<vector<float>> powerSpectrums;
    powerSpectrums.reserve(frameSplitSignal.size());

    pocketfft::shape_t  shape_in   = { (size_t)n_FFT };
    pocketfft::stride_t stride_in  = { sizeof(float) };
    pocketfft::stride_t stride_out = { sizeof(complex<float>) };

    vector<complex<float>> fftOutput(spectrumSize);

    vector<float> inputBuffer(n_FFT, 0.0f);
    cout << "frame size:" << frameSplitSignal[1].size() << endl;
    for (const auto& frame : frameSplitSignal) {
        
        size_t copyLen = min((size_t)n_FFT, frame.size());
        std::copy(frame.begin(), frame.begin() + copyLen, inputBuffer.begin());

        if (copyLen < n_FFT) {
            std::fill(inputBuffer.begin() + copyLen, inputBuffer.end(), 0.0f);
        }

        pocketfft::r2c<float>(
            shape_in, stride_in, stride_out, 0, true,
            inputBuffer.data(), fftOutput.data(), 1.0f
        );

        // --- 步骤 C: 计算功率谱 (Power Spectrum) ---
        // 公式: P = |X[k]|^2 / N
        vector<float> currentPowerSpec(spectrumSize);
        float nFloat = (float)n_FFT;

        for (size_t k = 0; k < spectrumSize; ++k) {
            float real = fftOutput[k].real();
            float imag = fftOutput[k].imag();
            

            currentPowerSpec[k] = (real * real + imag * imag) / nFloat;
        }
        // 存入结果矩阵
        powerSpectrums.push_back(std::move(currentPowerSpec));
    }

    
    // 此时 powerSpectrums 就是你需要的“功率谱”，接下来可以传给 Mel 滤波器组了
    // cout << "功率谱计算完成，总帧数: " << powerSpectrums.size() << endl;
    vector<vector<float>>melFilters = buildMELFilters(44100,n_melFilter);

    saveToBinary(powerSpectrums,"power.bin");

    vector<vector<float>>result(powerSpectrums.size());

    cout << "powerSpectrums size:" << powerSpectrums.size() << " " << powerSpectrums[1].size() << endl;
    cout << "melFilters size:" << melFilters.size() << " " << melFilters[1].size() << endl;

    for(int i = 0; i < powerSpectrums.size(); ++i)
    {
        for(int j = 0; j < n_melFilter; ++j)
        {
            float energy = 0;
            for(int k = 0; k < 1025; ++k)
            {
                energy += powerSpectrums[i][k] * melFilters[j][k];
                // cout << val << " " << powerSpectrums[i].at(k) << " " << melFilters[j].at(k) << endl;
            }
            if(energy < 1e-10) energy = 1e-10;
            result[i].push_back(logf(energy));   
            
        }
    }
    cout << "frame:" << result.size() << " length:" << result.at(1).size() << endl;
    saveToBinary(result,"result.bin");

    getchar();
    return 0;
}