# snoreDetectAlgorithm 💤

基于K-Means聚类的鼾声识别，旨在识别输入音频片段中是否包含鼾声。

## 📂 项目结构

| 文件/目录 | 说明 |
| :--- | :--- |
| `buildMFCC.cpp` | 提取声音特征 |
| `detectSnore.cpp` | 检测代码 |
| `main.cpp` | 测试代码  |
| `detect_snore_vq.m` | Matlab推理测试脚本 |
| `readMFCC.m` | 读取保存的MFCC文件并计算K-Means参数 |
| `saveCodebook.m` | 保存聚类结果为h文件 |
| `vq_codebooks.h`| 示例聚类结果h文件 |



## 🚀 使用流程

## 1) MFCC集创建
* 使用函数void saveAllMFCC()获取MFCC结果的.bin文件
* 执行readMFCC.m读取MFCC结果的.bin文件获取C_snore与C_noise
* （可选）执行detect_snore_vq.m测试C_snore与C_noise的效果
---

## 2) 识别

* 使用函数bool detectSnore()获取鼾声结果
---
