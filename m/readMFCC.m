clear
clc

rawData_0 = {};
rawData_1 = {};

for label=0:1
    for i = 0:499
        filename = "./mfcc/" + string(label) + "/" + string(0) + "_" + i + ".bin";
        fid = fopen(filename, 'r');
        if fid == -1
            disp(filename)
            error('无法打开文件');
        end
        % 1. 读取头部信息
        rows = fread(fid, 1, 'int32');
        cols = fread(fid, 1, 'int32');
        raw_data = fread(fid, [cols, rows], 'float');
    
        % 3. 关闭文件
        fclose(fid);
        
        % 4. 转置回来，变成 [rows, cols] (即 [帧数, 特征数])
        mfcc_data = raw_data';
    
        if(label == 0) rawData_0{end+1} = mfcc_data; end
        if(label == 1) rawData_1{end+1} = mfcc_data; end
    
    end
end

mfccData_0 = zeros(65,13,449);
mfccData_1 = zeros(65,13,449);

for i = 1:449
    mfccData_0(:,:,i) = cell2mat(rawData_0(i));
    mfccData_1(:,:,i) = cell2mat(rawData_1(i));
end

rawData_0 = cell2mat(rawData_0');
rawData_1 = cell2mat(rawData_1');

% return;
K_snore = 32;
K_noise = 64; % 噪音更复杂，给更多中心点

% 训练鼾声码书
% C_snore 是 32x13 的矩阵，每一行是一个中心向量
[idx_s, C_snore] = kmeans(rawData_1, K_snore, 'MaxIter', 1000);

% 训练噪音码书
% C_noise 是 64x13 的矩阵
[idx_n, C_noise] = kmeans(rawData_0, K_noise, 'MaxIter', 1000);

% 保存这两个矩阵，这就是你的"模型"
save('vq_codebooks.mat', 'C_snore', 'C_noise');




