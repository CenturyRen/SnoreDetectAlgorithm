import librosa
import librosa.display
import numpy as np
import matplotlib.pyplot as plt
import os
def process_wav_file(file_path):
    print(f"--- 正在处理文件: {file_path} ---")
    
    # ---------------------------------------------------------
    # 1. 加载并重采样 (Load & Resample)
    # ---------------------------------------------------------
    # 原始采样率是 44100Hz，但为了匹配论文和模型 (32x32)，
    # 我们必须在加载时强制重采样为 16000Hz (sr=16000)。
    # 鼾声主要集中在低频，16k 足够且能减少计算量。
    target_sr = 16000
    y, sr = librosa.load(file_path, sr=target_sr, duration=1.0)
    
    # 确保数据长度正好是 1秒 (16000点)
    # 如果文件略短，进行填充；如果略长，进行裁剪
    target_length = int(target_sr * 1.0)
    if len(y) < target_length:
        y = np.pad(y, (0, target_length - len(y)))
    else:
        y = y[:target_length]
        
    # print(f"原始采样率: 44100Hz -> 重采样后: {sr}Hz")
    # print(f"数据形状: {y.shape} (应为 16000)")

    # ---------------------------------------------------------
    # 2. 设置 MFCC 参数 (对应 32x32 输出)
    # ---------------------------------------------------------
    n_mfcc = 32      # 高度 (Y轴特征数)
    n_fft = 512      # FFT 窗口大小
    win_length = int(0.030 * sr) # 30ms 帧长
    
    # 关键计算：为了让 16000 个点最终变成宽度 32
    # Hop Length = 总点数 / 目标宽度 = 16000 / 32 = 500
    hop_length = 500 

    # ---------------------------------------------------------
    # 3. 提取特征
    # ---------------------------------------------------------
    mfcc = librosa.feature.mfcc(
        y=y, 
        sr=sr, 
        n_mfcc=n_mfcc, 
        n_mels=32,      
        n_fft=n_fft, 
        win_length=win_length, 
        hop_length=hop_length,
        center=True,
        # --- 新增参数以匹配论文 ---
        fmin=0,       # 论文暗示从低频开始覆盖
        fmax=8000,    # 明确覆盖到 8kHz (16k采样率的一半)
        htk=True      # 【重要】开启 HTK 模式以匹配 "linear up to 1kHz" 的描述
    )
    
    # ---------------------------------------------------------
    # 4. 尺寸对齐 (Hard Check)
    # ---------------------------------------------------------
    # 确保输出严格是 32x32
    target_width = 32
    if mfcc.shape[1] > target_width:
        mfcc = mfcc[:, :target_width]
    elif mfcc.shape[1] < target_width:
        mfcc = np.pad(mfcc, ((0, 0), (0, target_width - mfcc.shape[1])))
        
    # print(f"MFCC 矩阵尺寸: {mfcc.shape}")

    # ---------------------------------------------------------
    # 5. 标准化 (Normalization) - 对模型训练至关重要
    # ---------------------------------------------------------
    # 简单的 Z-score 标准化
    if np.std(mfcc) != 0:
        mfcc = (mfcc - np.mean(mfcc)) / np.std(mfcc)
    
    # ---------------------------------------------------------
    # 6. 可视化检查
    # ---------------------------------------------------------
    # plt.figure(figsize=(5, 4))
    # librosa.display.specshow(mfcc, sr=sr, hop_length=hop_length, x_axis='time')
    # plt.colorbar(label='Normalized Value')
    # plt.title(f'MFCC (32x32) from {file_path}')
    # plt.tight_layout()
    # plt.show()

    # ---------------------------------------------------------
    # 7. 格式化为模型输入
    # ---------------------------------------------------------
    # 增加 Batch 和 Channel 维度 -> (1, 32, 32, 1)
    model_input = mfcc[np.newaxis, ..., np.newaxis]
    # print(f"最终喂给模型的 Tensor 形状: {model_input.shape}")
    
    return model_input


def create_dataset(dataset_path):
    X = [] # 存放特征
    y = [] # 存放标签
    
    # 类别映射：文件夹名 -> 标签数字
    categories = {'0': 0, '1': 1}

    for folder, label in categories.items():
        folder_path = os.path.join(dataset_path, folder)
        if not os.path.exists(folder_path):
            print(f"警告: 路径不存在 {folder_path}")
            continue
            
        print(f"正在处理类别: {folder} (Label: {label})...")
        
        files = [f for f in os.listdir(folder_path) if f.endswith('.wav')]
        
        for file in files:
            file_path = os.path.join(folder_path, file)
            
            # 调用你的提取函数
            features = process_wav_file(file_path)
            
            if features is not None:
                X.append(features)
                y.append(label)

    # 转换为 Numpy 数组
    X = np.array(X)
    y = np.array(y)
    
    print(f"\n处理完成!")
    print(f"特征矩阵 X 形状: {X.shape}") # 应该是 (N, 32, 32, 1)
    print(f"标签向量 y 形状: {y.shape}") # 应该是 (N,)
    
    return X, y


# 运行
# 请确保当前目录下有 0_0.wav 文件
try:
    
    X_data, y_data = create_dataset("SnoringDataset")
    save_path = "snore_dataset.npz"
    np.savez(save_path, X=X_data, y=y_data)
    print(f"数据已保存至: {save_path}")


except FileNotFoundError:
    print("错误找不到路径。")