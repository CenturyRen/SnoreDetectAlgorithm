#include<stdio.h>
#include<vector>
#include"buildMFCC.h"
#include"vq_codebooks.h"
bool detectSnore(const std::vector<std::vector<float>>& mfccFeatures,const std::vector<std::vector<float>>& C_snore, const std::vector<std::vector<float>>& C_noise);
bool detectSnore(const std::vector<std::vector<float>>& mfccFeatures);
bool detectSnore(AudioData originSignal);
bool detectSnore(const std::vector<float>& sound);