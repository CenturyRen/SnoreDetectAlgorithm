
#include <iostream>
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

// 计算两个向量之间的欧氏距离
float euclidean_dist(const vector<float>& a, const vector<float>& b)
{
    if(a.size() != b.size()) return -1.0f;
    float sum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sum;

}

std::vector<float> pdist2(const std::vector<float>&a, const std::vector<std::vector<float>>& B)
{
    size_t m = B.size();
    std::vector<float> D(m, 0.0f);

    for (size_t j = 0; j < m; ++j) {
        D[j] = euclidean_dist(a, B[j]);
    }
    return D;
}

float total_dist(const std::vector<std::vector<float>>& mfccFeatures, const std::vector<std::vector<float>>& C)
{
    float total_distance = 0.0f;
    for (const auto& feature : mfccFeatures) {
        std::vector<float> distances = pdist2(feature, C);
        // 找到最小距离
        float min_distance = *std::min_element(distances.begin(), distances.end());
        total_distance += min_distance;
    }
    return total_distance;

}


bool detectSnore(const vector<vector<float>>& mfccFeatures,const vector<vector<float>>& C_snore, const vector<vector<float>>& C_noise) 
{
    
    auto N   = mfccFeatures.size();
    auto dim = mfccFeatures[0].size();

    float total_dist_snore = total_dist(mfccFeatures, C_snore)/N;
    float total_dist_noise = total_dist(mfccFeatures, C_noise)/N;

    cout << "N:" << N << " dim:" << dim << endl;
    return total_dist_snore < total_dist_noise ? true : false;
}

bool detectSnore(const vector<vector<float>>& mfccFeatures)
{
    vector<vector<float>> C_snore = arrayToVector(snore_codebook); // 假设已经定义并填充了 C_snore  
    vector<vector<float>> C_noise = arrayToVector(noise_codebook); // 假设已经定义并填充了 C_noise
    return detectSnore(mfccFeatures,C_snore,C_noise);
}

bool detectSnore(AudioData originSignal)
{
    vector<vector<float>>result = buildMFCC(originSignal,0.030,0.015,32,2048,13);
    return detectSnore(result);
}
bool detectSnore(const std::vector<float>& sound)
{
    return detectSnore(buildMFCC(sound));
}