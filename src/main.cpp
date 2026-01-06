#define _USE_MATH_DEFINES
#define DR_WAV_IMPLEMENTATION 

#include<stdio.h>
#include<pocketfft_hdronly.h>
#include "dr_wav.h"
// #include"buildMFCC.h"
#include"detectSnore.h"
using namespace std;


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


/**
 * 创建MFCC集
 */
void saveAllMFCC()
{
    
    for(int label=0; label<=1; ++label)
    {
        std::string root = "D:/Files/PyPrj/SnoreDetectModel/SnoringDataset/" + std::to_string(label) + "/";
        for (int i = 0; i < 449; ++i) { // 第二位数字
            std::string filename = std::to_string(label) + "_" + std::to_string(i) + ".wav";
                
            std::string pathToPath = root + filename;

            AudioData originSignal = readWav(pathToPath.c_str());
            if (originSignal.soungAudio.empty()) {
                printf("error 1");
                return;
            }
            vector<vector<float>>result = buildMFCC(originSignal,0.030,0.015,32,2048,13);
            saveToBinary(result,std::to_string(label) + "_" + std::to_string(i) + ".bin");
            cout << "saved file:" << std::to_string(label) + "_" + std::to_string(i) + ".bin" << endl;
        }
    }
}

int main()
{

    // saveAllMFCC();
    // return 0;
    int count = 0;
    std::string root = "D:/Files/PyPrj/SnoreDetectModel/SnoringDataset/1/";
    for (int i = 450; i < 499; ++i) { // 第二位数字
            // 拼接文件名
        std::string filename = std::to_string(1) + "_" + std::to_string(i) + ".wav";
            
        std::string pathToPath = root + filename;

        AudioData originSignal = readWav(pathToPath.c_str());
        if (originSignal.soungAudio.empty()) {
            printf("error 1");
            return -1;
        }


        vector<float>momoSignal;
        if(originSignal.channels==2)
            momoSignal = stereoToMono(originSignal.soungAudio);
        else
            momoSignal = originSignal.soungAudio;

        auto isSnore = detectSnore(momoSignal);
        cout << "process file:" << filename << " " << isSnore << endl;
        if(!isSnore)
        {
            count++;
        } 
    }

    cout << "wrong count:" << count << endl;
    getchar();
    return 0;
}