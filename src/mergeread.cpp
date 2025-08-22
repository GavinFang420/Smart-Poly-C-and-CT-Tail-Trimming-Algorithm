#include "smarttrim.h"
#include "mergeread.h"
#include <algorithm>
#include <iostream>
#include <cmath>

// 修正的动态递进检测算法 - 只分析merge后的C tail
TrimResult SmartTrimmer::findOptimalTrimPositions(
    const std::string& r1_seq, 
    const std::string& r2_seq
) {
    TrimResult result;
    
    if (r1_seq.empty() || r2_seq.empty()) {
        return result;
    }
    
    // 先尝试merge
    ReadMerger merger(10, 3, 0.3);
    MergeResult merge_result = merger.mergeReads(r1_seq, r2_seq);
    
    if (merge_result.merged) {
        // Merge成功，分析merged序列尾部的polyC tail
        return analyzeMergedSequence(merge_result, r1_seq, r2_seq);
    } else {
        // Merge失败，直接分析R2头部的polyG tail (作为fallback)
        return analyzeR2PolyGTail(r2_seq);
    }
}

TrimResult SmartTrimmer::analyzeMergedSequence(
    const MergeResult& merge_result,
    const std::string& original_r1_seq,
    const std::string& original_r2_seq
) {
    TrimResult result;
    std::string merged_seq = merge_result.sequence;
    
    // 从merged序列末尾开始动态检测target tail (支持C, CT, CG等)
    double current_score = params.initial_score;
    int best_cutoff = 0;
    double best_cutoff_score = params.initial_score;
    double cutoff_next_score = 0.0;
    bool has_next_score = false;
    bool in_death_mode = false;
    int death_mode_start = -1;
    int consecutive_good_bases = 0;
    int consecutive_target_count = 0;  // 连续目标碱基计数
    
    int max_check_length = std::min(params.window_size, (int)merged_seq.length());
    
    for (int i = 0; i < max_check_length; i++) {
        int pos = merged_seq.length() - 1 - i; // 从末尾开始
        char base = std::toupper(merged_seq[pos]);
        
        // 使用新的权重计算系统
        double distance_weight = params.calculateDistanceWeight(i);
        
        // 连续目标碱基权重计算
        double consecutive_weight = 1.0;
        if (params.isTargetBase(base)) {
            consecutive_target_count++;
            consecutive_weight = params.calculateConsecutiveWeight(consecutive_target_count);
            consecutive_good_bases++;
        } else {
            consecutive_target_count = 0; // 重置连续计数
            consecutive_good_bases = 0;
        }
        
        // 最终权重 = 距离权重 × 连续权重
        double final_weight = distance_weight * consecutive_weight;
        
        if (params.isTargetBase(base)) {
            // 目标碱基加分
            double base_score = params.getBaseScore(base);
            current_score += base_score * final_weight;
            
            // 濒死模式检查复活
            if (in_death_mode && consecutive_good_bases >= 3) {
                in_death_mode = false;
                death_mode_start = -1;
                consecutive_good_bases = 0;
            }
        } else {
            // 非目标碱基扣分
            double base_score = params.getBaseScore(base);
            double penalty_multiplier = 1.0 + 2.0 * distance_weight; // 1.0到3.0之间
            current_score += base_score * penalty_multiplier;
            
            // 记录cutoff后第一个非目标碱基的分数
            if (best_cutoff > 0 && !has_next_score) {
                cutoff_next_score = current_score;
                has_next_score = true;
            }
            
            // 检查是否进入濒死模式
            if (current_score < 0.0 && !in_death_mode) {
                in_death_mode = true;
                death_mode_start = i;
                consecutive_target_count = 0; // 重置连续计数，重新开始累积
            }
        }
        
        // 如果分数还是正的，更新最佳cutoff
        if (current_score >= 0.0) {
            best_cutoff = i + 1;
            best_cutoff_score = current_score;
        }
        
        // 濒死模式超时检查 - 如果连续5个非目标碱基就坠机
        if (in_death_mode && (i - death_mode_start) > 5) {
            break; // 坠机，停止检测
        }
    }
    
    if (best_cutoff > 0) {
        // 映射到原始reads的trim位置
        std::pair<int, int> trim_positions = mapMergedPositionToOriginal(
            best_cutoff, merge_result, original_r1_seq, original_r2_seq
        );
        
        result.r1_trim_pos = trim_positions.first;
        result.r2_trim_pos = trim_positions.second;
        result.final_score = best_cutoff_score;
        
        // 生成score详情
        if (has_next_score) {
            result.score_detail = std::to_string(best_cutoff_score) + " -> " + std::to_string(cutoff_next_score);
        } else {
            result.score_detail = std::to_string(best_cutoff_score) + " -> (no next)";
        }
        
        result.is_valid = true;
    }
    
    return result;
}

TrimResult SmartTrimmer::analyzeR2PolyGTail(const std::string& r2_seq) {
    TrimResult result;
    
    // Fallback: 当merge失败时，分析R2头部的target tail
    double current_score = params.initial_score;
    int best_cutoff = 0;
    double best_cutoff_score = params.initial_score;
    double cutoff_next_score = 0.0;
    bool has_next_score = false;
    bool in_death_mode = false;
    int death_mode_start = -1;
    int consecutive_good_bases = 0;
    int consecutive_target_count = 0;  // 连续目标碱基计数
    
    int max_check_length = std::min(params.window_size, (int)r2_seq.length());
    
    for (int i = 0; i < max_check_length; i++) {
        char base = std::toupper(r2_seq[i]);
        
        // 距离衰减：从R2头部开始
        double distance_weight = params.calculateDistanceWeight(i);
        
        // 连续目标碱基衰减
        double consecutive_weight = 1.0;
        if (params.isTargetBase(base)) {
            consecutive_target_count++;
            consecutive_weight = params.calculateConsecutiveWeight(consecutive_target_count);
            consecutive_good_bases++;
        } else {
            consecutive_target_count = 0;
            consecutive_good_bases = 0;
        }
        
        double final_weight = distance_weight * consecutive_weight;
        
        if (params.isTargetBase(base)) {
            double base_score = params.getBaseScore(base);
            current_score += base_score * final_weight;
            
            // 濒死模式检查复活
            if (in_death_mode && consecutive_good_bases >= 3) {
                in_death_mode = false;
                death_mode_start = -1;
                consecutive_good_bases = 0;
            }
        } else {
            // 非目标碱基扣分
            double base_score = params.getBaseScore(base);
            double penalty_multiplier = 1.0 + 2.0 * distance_weight;
            current_score += base_score * penalty_multiplier;
            
            if (best_cutoff > 0 && !has_next_score) {
                cutoff_next_score = current_score;
                has_next_score = true;
            }
            
            // 检查是否进入濒死模式
            if (current_score < 0.0 && !in_death_mode) {
                in_death_mode = true;
                death_mode_start = i;
                consecutive_target_count = 0; // 重置连续计数
            }
        }
        
        // 如果分数还是正的，更新最佳cutoff
        if (current_score >= 0.0) {
            best_cutoff = i + 1;
            best_cutoff_score = current_score;
        }
        
        // 濒死模式超时检查
        if (in_death_mode && (i - death_mode_start) > 5) {
            break; // 坠机
        }
    }
    
    if (best_cutoff > 0) {
        result.r1_trim_pos = 0; // 不trim R1
        result.r2_trim_pos = best_cutoff; // 从R2开头trim
        result.final_score = best_cutoff_score;
        
        if (has_next_score) {
            result.score_detail = std::to_string(best_cutoff_score) + " -> " + std::to_string(cutoff_next_score);
        } else {
            result.score_detail = std::to_string(best_cutoff_score) + " -> (no next)";
        }
        
        result.is_valid = true;
    }
    
    return result;
}

// 修正的映射函数 - 考虑生物学对应关系
std::pair<int, int> SmartTrimmer::mapMergedPositionToOriginal(
    int cut_length_from_end, 
    const MergeResult& merge_result,
    const std::string& original_r1_seq,
    const std::string& original_r2_seq
) {
    std::string merged_seq = merge_result.sequence;
    int merged_length = merged_seq.length();
    
    // 获取要trim的polyC region
    std::string polyc_region = merged_seq.substr(merged_length - cut_length_from_end, cut_length_from_end);
    
    int r1_trim = 0;
    int r2_trim = 0;
    
    // 关键：分析polyC tail的生物学来源
    // 检查原始R1末尾的polyC含量
    int r1_end_c_count = 0;
    int check_length_r1 = std::min(cut_length_from_end, (int)original_r1_seq.length());
    if (check_length_r1 > 0) {
        std::string r1_end = original_r1_seq.substr(original_r1_seq.length() - check_length_r1);
        for (char base : r1_end) {
            if (std::toupper(base) == 'C' || std::toupper(base) == 'N') {
                r1_end_c_count++;
            }
        }
    }
    
    // 检查原始R2开头的polyG含量（对应polyC）
    int r2_start_g_count = 0;
    int check_length_r2 = std::min(cut_length_from_end, (int)original_r2_seq.length());
    if (check_length_r2 > 0) {
        std::string r2_start = original_r2_seq.substr(0, check_length_r2);
        for (char base : r2_start) {
            if (std::toupper(base) == 'G' || std::toupper(base) == 'N') {
                r2_start_g_count++;
            }
        }
    }
    
    // 生物学对应关系：
    // 如果R1末尾有polyC且R2开头有polyG，它们在merge后都变成同一段polyC tail
    // 我们需要同时trim它们！
    
    if (r1_end_c_count >= 3 && r2_start_g_count >= 3) {
        // 两者都有明显的poly特征 - 同时trim
        r1_trim = r1_end_c_count;  // trim R1末尾的所有C
        r2_trim = r2_start_g_count; // trim R2开头的所有G
    } else if (r1_end_c_count >= 3) {
        // 只有R1有polyC - 只trim R1
        r1_trim = r1_end_c_count;
        r2_trim = 0;
    } else if (r2_start_g_count >= 3) {
        // 只有R2有polyG - 只trim R2
        r1_trim = 0;
        r2_trim = r2_start_g_count;
    } else {
        // 都没有明显特征 - 保守策略
        // 但如果算法检测到了polyC tail，说明确实有问题
        // 按检测到的长度分配
        if (r1_end_c_count > 0 && r2_start_g_count > 0) {
            r1_trim = r1_end_c_count;
            r2_trim = r2_start_g_count;
        } else if (r1_end_c_count > 0) {
            r1_trim = cut_length_from_end;
            r2_trim = 0;
        } else {
            r1_trim = 0;
            r2_trim = cut_length_from_end;
        }
    }
    
    // 确保不超出原始序列长度
    r1_trim = std::max(0, std::min(r1_trim, (int)original_r1_seq.length()));
    r2_trim = std::max(0, std::min(r2_trim, (int)original_r2_seq.length()));
    
    return std::make_pair(r1_trim, r2_trim);
}

// trim函数
std::pair<std::string, std::string> SmartTrimmer::trimReads(
    const std::string& r1_seq,
    const std::string& r2_seq,
    const TrimResult& result
) {
    if (!result.is_valid) {
        return std::make_pair(r1_seq, r2_seq);
    }
    
    std::string trimmed_r1 = r1_seq;
    std::string trimmed_r2 = r2_seq;
    
    // 从R1末尾trim
    if (result.r1_trim_pos > 0 && result.r1_trim_pos < r1_seq.length()) {
        trimmed_r1 = r1_seq.substr(0, r1_seq.length() - result.r1_trim_pos);
    }
    
    // 从R2开头trim
    if (result.r2_trim_pos > 0 && result.r2_trim_pos < r2_seq.length()) {
        trimmed_r2 = r2_seq.substr(result.r2_trim_pos);
    }
    
    return std::make_pair(trimmed_r1, trimmed_r2);
}

// 参数矩阵生成 - 支持CT tail和PolyC tail测试
std::vector<TrimParams> SmartTrimmer::generateParameterMatrix() {
    std::vector<TrimParams> param_list;
    
    std::vector<int> window_sizes = {15, 20, 25, 30};
    std::vector<double> initial_scores = {0.0, -5.0, -10.0};
    
    // PolyC配置矩阵
    std::vector<double> c_scores = {8.0, 10.0, 12.0, 15.0};
    std::vector<double> penalty_scores = {-15.0, -20.0, -25.0};
    
    for (int ws : window_sizes) {
        for (double c_score : c_scores) {
            for (double penalty : penalty_scores) {
                for (double init_score : initial_scores) {
                    TrimParams params = TrimParams::createPolyCConfig(ws);
                    params.c_score = c_score;
                    params.a_score = penalty;
                    params.g_score = penalty;
                    params.t_score = penalty;
                    params.initial_score = init_score;
                    param_list.push_back(params);
                }
            }
        }
    }
    
    // CT tail配置矩阵
    std::vector<double> t_scores = {6.0, 8.0, 10.0}; // T分数通常比C分数低
    std::vector<double> ct_ratios = {0.6, 0.8, 1.0}; // T/C 比例
    
    for (int ws : window_sizes) {
        for (double c_score : c_scores) {
            for (double t_ratio : ct_ratios) {
                for (double penalty : penalty_scores) {
                    for (double init_score : initial_scores) {
                        TrimParams params = TrimParams::createCTConfig(ws);
                        params.c_score = c_score;
                        params.t_score = c_score * t_ratio;
                        params.a_score = penalty;
                        params.g_score = penalty;
                        params.initial_score = init_score;
                        param_list.push_back(params);
                    }
                }
            }
        }
    }
    
    return param_list;
}
