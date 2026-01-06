#include"buildMFCC.h"

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

// 生成 DCT 变换矩阵
// nFilters: 滤波器个数 (32)
// nCepstral: 保留的 MFCC 系数个数 (通常取 13)
std::vector<std::vector<float>> createDCTMatrix(int nFilters, int nCepstral) {
    std::vector<std::vector<float>> dctMat(nCepstral, std::vector<float>(nFilters));
    
    // 正交归一化系数 (Orthogonal Normalization)
    // 这能保证变换后的数值范围不会爆炸
    float normalizer = sqrt(2.0f / nFilters);
    
    for (int k = 0; k < nCepstral; ++k) {
        for (int n = 0; n < nFilters; ++n) {
            // DCT-II 公式
            float angle = M_PI / nFilters * (n + 0.5f) * k;
            float value = cos(angle) * normalizer;
            
            // 修正 k=0 时的归一化系数 (1/sqrt(2))
            if (k == 0) {
                value *= sqrt(0.5f);
            }
            
            dctMat[k][n] = value;
        }
    }
    return dctMat;
}

// 执行 DCT 变换
// 输入: Log Fbank 能量 [32]
// 输出: MFCC 系数 [13]
std::vector<float> applyDCT(const std::vector<float>& logEnergies, const std::vector<std::vector<float>>& dctMat) {
    int nCepstral = dctMat.size();
    int nFilters = dctMat[0].size();
    
    std::vector<float> mfcc(nCepstral, 0.0f);
    
    for (int k = 0; k < nCepstral; ++k) {
        float sum = 0.0f;
        for (int n = 0; n < nFilters; ++n) {
            sum += logEnergies[n] * dctMat[k][n];
        }
        mfcc[k] = sum;
    }
    return mfcc;
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

// vector<vector<float>>buildMFFCC(const <vector<float>& sound, int sampleRate = 44100,float stepTime = 0.030 ,float stepDuration = 0.015, int n_melFilter = 32, int n_FFT = 2048, int numCepstral = 13)
vector<vector<float>>buildMFCC(const vector<float>& sound, int sampleRate ,float stepTime,float stepDuration , int n_melFilter, int n_FFT, int numCepstral)

{

    int win_Length = int(stepTime * sampleRate);// 窗口宽度（30ms的宽度）

    int stepLength = int(stepDuration * sampleRate);
    int n_MFCC = (sound.size() - win_Length) / stepLength + 1;
    // 预加强
    vector<float> preE_Signal = preEmphasis(sound);

    // 分帧
    vector<vector<float>> frameSplitSignal;
    frameSplitSignal.reserve(n_MFCC); // 预留外层空间，提高效率

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

    size_t spectrumSize = n_FFT / 2 + 1;
    vector<vector<float>> powerSpectrums;
    powerSpectrums.reserve(frameSplitSignal.size());

    pocketfft::shape_t  shape_in   = { (size_t)n_FFT };
    pocketfft::stride_t stride_in  = { sizeof(float) };
    pocketfft::stride_t stride_out = { sizeof(complex<float>) };

    vector<complex<float>> fftOutput(spectrumSize);

    vector<float> inputBuffer(n_FFT, 0.0f);
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
    vector<vector<float>>result(powerSpectrums.size());

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

    // 1. 创建 DCT 矩阵 (只做一次，放循环外)
    vector<vector<float>> dctMatrix = createDCTMatrix(32, numCepstral);
    // 2. 准备最终容器
    vector<vector<float>> mfccResult;
    mfccResult.reserve(result.size());
    
    for (const auto& fbankFrame : result) {
        vector<float> mfccFrame = applyDCT(fbankFrame, dctMatrix);
        mfccResult.push_back(mfccFrame);
    }

    return mfccResult;
}
/**创建mfcc
 * originSignal:读取出的音频结构体
 * stepTime:窗口长度（ms）
 * stepDuration:重叠长度（ms）
 * n_melFilter:mfl滤波器数量
 * n_FFT:FFT滤波器长度
 * numCepstral:要保留的特征数量
 */
vector<vector<float>>buildMFCC(AudioData &originSignal,float stepTime,float stepDuration,int n_melFilter,int n_FFT,int numCepstral)
{
    vector<float>momoSignal;
    // 双声道转单声道
    if(originSignal.channels==2)
        momoSignal = stereoToMono(originSignal.soungAudio);
    else
        momoSignal = originSignal.soungAudio;
    // 定义MFCC信息
    unsigned int sampleRate = originSignal.sampleRate;

    int win_Length = int(stepTime * sampleRate);// 窗口宽度（30ms的宽度）

    int stepLength = int(stepDuration * sampleRate);
    int n_MFCC = (momoSignal.size() - win_Length) / stepLength + 1;
    // 预加强
    vector<float> preE_Signal = preEmphasis(momoSignal);

    // 分帧
    vector<vector<float>> frameSplitSignal;
    frameSplitSignal.reserve(n_MFCC); // 预留外层空间，提高效率

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
    vector<vector<float>>result(powerSpectrums.size());

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

    // 1. 创建 DCT 矩阵 (只做一次，放循环外)
    vector<vector<float>> dctMatrix = createDCTMatrix(32, numCepstral);
    // 2. 准备最终容器
    vector<vector<float>> mfccResult;
    mfccResult.reserve(result.size());
    
    for (const auto& fbankFrame : result) {
        vector<float> mfccFrame = applyDCT(fbankFrame, dctMatrix);
        mfccResult.push_back(mfccFrame);
    }

    return mfccResult;

}


