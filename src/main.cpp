#define _USE_MATH_DEFINES
#define DR_WAV_IMPLEMENTATION 

#include<stdio.h>
#include<pocketfft_hdronly.h>
#include "dr_wav.h"
#include"buildMFCC.h"
#include"detectSnore.h"
using namespace std;


template <int R, int C>
std::vector<std::vector<float>> arrayToVector(const float (&arr)[R][C]) {
    std::vector<std::vector<float>> result;
    result.reserve(R);
    for (int i = 0; i < R; ++i) {
        result.emplace_back(arr[i], arr[i] + C);
    }
    return result;
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
    // std::cout << "save seccess:" << filename << std::endl;
}


int main()
{

    std::string root = "D:/Files/PyPrj/SnoreDetectModel/SnoringDataset/0/";
    for (int i = 450; i < 499; ++i) { // 第二位数字
            // 拼接文件名
        std::string filename = std::to_string(0) + "_" + std::to_string(i) + ".wav";
            
        std::string pathToPath = root + filename;

        AudioData originSignal = readWav(pathToPath.c_str());
        if (originSignal.soungAudio.empty()) {
            printf("error 1");
            return -1;
        }
        vector<vector<float>>result = buildMFCC(originSignal,0.030,0.015,32,2048,13);
        vector<vector<float>> C_snore = arrayToVector(snore_codebook); // 假设已经定义并填充了 C_snore  
        vector<vector<float>> C_noise = arrayToVector(noise_codebook); // 假设已经定义并填充了 C_noise

        auto isSnore = detectSnore(result, C_snore, C_noise);
        cout << "process file:" << filename << " " << isSnore << endl;
    }
    
    getchar();
    return 0;
}