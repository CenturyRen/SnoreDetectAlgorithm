function isSnore = detect_snore_vq(test_mfcc, C_snore, C_noise)
    % test_mfcc: 13 x N 的矩阵
    % C_snore: K_snore x 13
    % C_noise: K_noise x 13
    
    % frames = test_mfcc; % 转置为 N x 13
    [N, dim] = size(test_mfcc);
    
    total_dist_snore = 0;
    total_dist_noise = 0;
    
    % 对每一帧计算失真
    for t = 1:N
        current_frame = test_mfcc(t, :);
        
        % 1. 找离鼾声码书最近的距离
        % pdist2 计算欧氏距离
        d_s = min(pdist2(current_frame, C_snore)); 
        total_dist_snore = total_dist_snore + d_s;
        
        % 2. 找离噪音码书最近的距离
        d_n = min(pdist2(current_frame, C_noise));
        total_dist_noise = total_dist_noise + d_n;
    end
    
    % 计算平均失真
    avg_dist_snore = total_dist_snore / N;
    avg_dist_noise = total_dist_noise / N;
    
    % 判决逻辑
    % 可以在这里加一个加权系数 alpha，如果噪音漏检率高，可以适当降低噪音的距离权重
    if avg_dist_snore < avg_dist_noise
        isSnore = true;
    else
        isSnore = false;
    end
end