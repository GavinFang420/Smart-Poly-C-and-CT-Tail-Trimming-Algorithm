#include "smarttrim.h"
#include "mergeread.h"
#include <algorithm>
#include <iostream>
#include <cmath>

double SmartTrimmer::calculateProgressiveScore(const std::string& sequence, int start_pos, int end_pos) {
    if (start_pos >= end_pos || start_pos < 0 || end_pos > sequence.length()) {
        return params.initial_score;
    }
    
    double score = params.initial_score;
    int region_len = end_pos - start_pos;
    
    // 分析指定区域，使用Distance-Decay权重
    for (int i = 0; i < region_len; i++) {
        char base = sequence[start_pos + i];
        
        // 位置权重：Distance-Decay，越靠近开头权重越大（因为是poly-C tail）
        double position_weight = (i < params.position_weights.size()) ? 
            params.position_weights[region_len - 1 - i] : 1.0;  // 反向权重
        
        if (base == 'C' || base == 'c') {
            score += params.c_score * position_weight;
        } else {
            // A, G, T get penalty
            score += params.penalty_score * position_weight;
        }
    }
    
    return score;
}

TrimResult SmartTrimmer::findOptimalTrimPositions(
    const std::string& r1_seq, 
    const std::string& r2_seq
) {
    TrimResult result;
    
    if (r1_seq.empty() || r2_seq.empty()) {
        return result;  // Invalid input
    }
    
    // 直接分析R2开头的poly-C tail
    // 不需要merge，直接检测R2序列开头
    
    double best_score = params.initial_score;
    int best_cut_pos = 0;
    
    // 检测R2开头最多window_size长度的poly-C
    int max_check_length = std::min(params.window_size, (int)r2_seq.length());
    
    // 尝试不同的cut位置（从R2开头开始）
    for (int cut_pos = 1; cut_pos <= max_check_length; cut_pos++) {
            // 首先检查C含量是否达标
        std::string candidate_region = r2_seq.substr(0, cut_pos);
        int c_count = 0;
        for (char base : candidate_region) {
            if (base == 'C' || base == 'c') c_count++;
        }
        double c_ratio = (double)c_count / cut_pos;
        
        // C含量必须≥80%才考虑
        if (c_ratio < 0.8) continue;  // ← 添加这个检查
        // ... 剩余逻辑
        // 分析R2开头cut_pos长度的区域
        double score = calculateProgressiveScore(r2_seq, 0, cut_pos);
        
        if (score > best_score) {
            best_score = score;
            best_cut_pos = cut_pos;
        }
    }
    
    // 应用阈值检查
    if (best_score > 0.0 && best_cut_pos > 0) {
        result.r1_trim_pos = 0;  // 不trim R1
        result.r2_trim_pos = best_cut_pos;  // 从R2开头trim
        result.final_score = best_score;
        result.is_valid = true;
    }
    
    return result;
}

// 保留legacy函数以防需要
double SmartTrimmer::calculateProgressiveScore(const std::string& merged_seq) {
    return calculateProgressiveScore(merged_seq, 0, merged_seq.length());
}

std::pair<int, int> SmartTrimmer::mapToOriginalPositions(int merge_pos, int r1_len, int r2_len) {
    // Legacy function - 保持兼容性
    return std::make_pair(0, merge_pos);
}

TrimResult SmartTrimmer::findOptimalTrimPositions_Window(
    const std::string& r1_seq, 
    const std::string& r2_seq
) {
    // 这个函数现在和主函数逻辑相同，因为我们总是分析R2开头
    return findOptimalTrimPositions(r1_seq, r2_seq);
}

std::pair<int, int> SmartTrimmer::mapMergedPositionToOriginal(
    int cut_length_in_r2, 
    const MergeResult& merge_result,
    int original_r1_length,
    int original_r2_length
) {
    // 简化：直接返回R2的trim位置
    int r1_trim = 0;
    int r2_trim = std::min(cut_length_in_r2, original_r2_length);
    
    return std::make_pair(r1_trim, r2_trim);
}

std::pair<std::string, std::string> SmartTrimmer::trimReads(
    const std::string& r1_seq,
    const std::string& r2_seq,
    const TrimResult& result
) {
    if (!result.is_valid) {
        return std::make_pair(r1_seq, r2_seq);  // No trimming
    }
    
    // R1保持不变
    std::string trimmed_r1 = r1_seq;
    
    // 从R2开头trim指定长度
    std::string trimmed_r2 = r2_seq;
    if (result.r2_trim_pos > 0 && result.r2_trim_pos < r2_seq.length()) {
        trimmed_r2 = r2_seq.substr(result.r2_trim_pos);
    }
    
    return std::make_pair(trimmed_r1, trimmed_r2);
}

std::vector<TrimParams> SmartTrimmer::generateParameterMatrix() {
    std::vector<TrimParams> param_list;
    
    // Test different window sizes for R2 head analysis
    std::vector<int> window_sizes = {15, 20, 25, 30};
    
    // Test different penalty scores
    std::vector<double> penalty_scores = {-2.0, -3.0, -4.0, -5.0};
    
    // Test different initial scores
    std::vector<double> initial_scores = {0.0, -10.0, -20.0};
    
    for (int ws : window_sizes) {
        for (double penalty : penalty_scores) {
            for (double init_score : initial_scores) {
                TrimParams params(ws, init_score);
                params.penalty_score = penalty;
                param_list.push_back(params);
            }
        }
    }
    
    // Generate different weight decay functions for R2 head analysis
    for (int ws : window_sizes) {
        // Exponential decay (强调开头)
        TrimParams exp_params(ws, 0.0);
        for (int i = 0; i < ws; i++) {
            // 越靠近开头权重越大
            exp_params.position_weights[i] = std::exp(2.0 * (ws - 1 - i) / (ws - 1));
        }
        param_list.push_back(exp_params);
        
        // Quadratic decay (强调开头)
        TrimParams quad_params(ws, 0.0);
        for (int i = 0; i < ws; i++) {
            double x = (double)(ws - 1 - i) / (ws - 1);  // 反向
            quad_params.position_weights[i] = 1.0 + 2.0 * x * x;
        }
        param_list.push_back(quad_params);
    }
    
    std::cout << "Generated " << param_list.size() << " parameter configurations for R2 head poly-C analysis." << std::endl;
    return param_list;
}