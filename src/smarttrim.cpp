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

// 预检查函数 - 检测target base富集情况
double SmartTrimmer::performPreCheck(const std::string& sequence, bool is_tail_check) const {
    int check_length = std::min(25, (int)sequence.length());
    if (check_length == 0) return 0.0;
    
    int target_count = 0;
    
    if (is_tail_check) {
        // 检查尾部
        for (int i = 0; i < check_length; i++) {
            int pos = sequence.length() - 1 - i;
            char base = std::toupper(sequence[pos]);
            if (params.isTargetBase(base)) {
                target_count++;
            }
        }
    } else {
        // 检查头部
        for (int i = 0; i < check_length; i++) {
            char base = std::toupper(sequence[i]);
            if (params.isTargetBase(base)) {
                target_count++;
            }
        }
    }
    
    return (double)target_count / check_length;
}

// 强制trim条件检查
int SmartTrimmer::checkForceTrimmingCondition(const std::string& sequence, 
                                              double target_ratio, bool is_tail_check) const {
    if (target_ratio < params.force_trim_threshold) {
        return 0;
    }
    
    int check_length = std::min(params.max_force_trim_length, (int)sequence.length());
    
    if (is_tail_check) {
        // 从尾部开始找第一个非target base
        for (int i = 1; i <= check_length; i++) {
            int pos = sequence.length() - i;
            char base = std::toupper(sequence[pos]);
            if (!params.isTargetBase(base)) {
                return std::max(1, i - 1);  // 保留一个非target base
            }
        }
        return std::min(check_length, 15);  // 最多强制trim 15bp
    } else {
        // 从头部开始找第一个非target base
        for (int i = 1; i <= check_length; i++) {
            char base = std::toupper(sequence[i-1]);
            if (!params.isTargetBase(base)) {
                return std::max(1, i - 1);
            }
        }
        return std::min(check_length, 15);
    }
}

TrimResult SmartTrimmer::analyzeMergedSequence(
    const MergeResult& merge_result,
    const std::string& original_r1_seq,
    const std::string& original_r2_seq
) {
    TrimResult result;
    std::string merged_seq = merge_result.sequence;
    
    // 实时权重和概率计算
    double current_score = params.initial_score;
    int best_cutoff = 0;
    double best_cutoff_score = params.initial_score;
    double cutoff_next_score = 0.0;
    bool has_next_score = false;
    
    // 添加濒死状态控制
    bool in_death_mode = false;
    int death_mode_start = -1;
    int death_count = 0;
    const int max_death_count = 2;  // 濒死2次后停止
    
    int consecutive_ct_count = 0;
    int consecutive_ag_count = 0;
    double current_ct_weight = 1.0;
    double current_ag_penalty_multiplier = 1.0;
    
    int max_check_length = std::min(std::max(params.window_size, 50), (int)merged_seq.length());
    
    // 主循环：分数没cutoff就一直走
    for (int i = 0; i < max_check_length; i++) {
        int pos = merged_seq.length() - 1 - i; // 从末尾开始
        char base = std::toupper(merged_seq[pos]);
        
        // 距离权重
        double distance_weight = params.calculateDistanceWeight(i);
        
        // 判断碱基类型并更新权重
        bool is_ct = (base == 'C' || base == 'T');
        bool is_ag = (base == 'A' || base == 'G');
        
        if (is_ct) {
            // CT碱基：累加连续计数，应用衰减权重
            consecutive_ct_count++;
            consecutive_ag_count = 0; // 重置AG计数
            
            // CT权重衰减：每连续出现一次衰减
            current_ct_weight *= params.consecutive_ct_decay_rate;
            current_ag_penalty_multiplier = 1.0; // 重置AG惩罚倍数
            
            // CT碱基加分
            double base_score = params.getBaseScore(base);
            double final_weight = distance_weight * current_ct_weight;
            current_score += base_score * final_weight;
            
            // 濒死模式复活检查：连续3个CT可以复活
            if (in_death_mode && consecutive_ct_count >= 3) {
                in_death_mode = false;
                death_mode_start = -1;
            }
            
        } else if (is_ag) {
            // AG碱基：高惩罚且权重递增
            consecutive_ag_count++;
            consecutive_ct_count = 0; // 重置CT计数
            
            // AG惩罚权重递增：每次增加20%
            current_ag_penalty_multiplier *= (1.0 + params.ag_penalty_increment);
            current_ct_weight = 1.0; // 重置CT权重
            
            // AG碱基扣分 - 前8个位置扣分减半
            double base_penalty = params.ag_penalty_base;
            double final_penalty = distance_weight * current_ag_penalty_multiplier;
            
            // 前8个位置扣分减半
            if (i < 8) {
                final_penalty *= 0.5;
            }
            
            current_score += base_penalty * final_penalty;
            
            // 检查是否进入濒死模式
            if (current_score < -20.0 && !in_death_mode) {  // 濒死阈值
                in_death_mode = true;
                death_mode_start = i;
            }
            
        } else {
            // N或其他碱基
            consecutive_ct_count = 0;
            consecutive_ag_count = 0;
            current_ct_weight = 1.0;
            current_ag_penalty_multiplier = 1.0;
            
            double base_score = params.getBaseScore(base);
            current_score += base_score * distance_weight;
        }
        
        // 记录下一个分数
        if (best_cutoff > 0 && !has_next_score) {
            cutoff_next_score = current_score;
            has_next_score = true;
        }
        
        // 无条件更新最佳cutoff（分数没cutoff就一直走）
        int potential_cutoff = i + 1;
        
        // 检查概率和边界条件
        if (potential_cutoff > 0) {
            // 计算要切除的序列概率
            int cut_start = merged_seq.length() - potential_cutoff;
            int cut_end = merged_seq.length();
            double cut_probability = params.calculateCTTailProbability(merged_seq, cut_start, cut_end);
            
            // 检查切除边界
            bool valid_boundary = params.isValidCutBoundary(merged_seq, potential_cutoff, true);
            
            // 如果概率满足条件且边界有效，更新cutoff
            if (cut_probability >= params.min_ct_tail_probability && valid_boundary) {
                best_cutoff = potential_cutoff;
                best_cutoff_score = current_score;
            }
        }
        
        // 濒死模式超时检查
        if (in_death_mode && (i - death_mode_start) > 10) {  // 10次后超时
            death_count++;
            if (death_count >= max_death_count) {
                break; // 濒死2次后彻底停止
            }
            // 重置濒死状态，继续下一轮濒死
            in_death_mode = false;
            death_mode_start = -1;
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
        
        // 生成详细信息
        int cut_start = merged_seq.length() - best_cutoff;
        int cut_end = merged_seq.length();
        double final_probability = params.calculateCTTailProbability(merged_seq, cut_start, cut_end);
        
        result.score_detail = "score=" + std::to_string(best_cutoff_score) + 
                             ",prob=" + std::to_string(final_probability) +
                             ",len=" + std::to_string(best_cutoff);
        
        result.is_valid = true;
    }
    
    return result;
}

TrimResult SmartTrimmer::analyzeR2PolyGTail(const std::string& r2_seq) {
    TrimResult result;
    
    // R2头部实时权重和概率计算
    double current_score = params.initial_score;
    int best_cutoff = 0;
    double best_cutoff_score = params.initial_score;
    double cutoff_next_score = 0.0;
    bool has_next_score = false;
    
    // 添加濒死状态控制
    bool in_death_mode = false;
    int death_mode_start = -1;
    int death_count = 0;
    const int max_death_count = 2;  // 濒死2次后停止
    
    int consecutive_ga_count = 0;
    int consecutive_ct_count = 0;
    double current_ga_weight = 1.0;
    double current_ct_penalty_multiplier = 1.0;
    
    int max_check_length = std::min(std::max(params.window_size, 50), (int)r2_seq.length());
    
    // 主循环：分数没cutoff就一直走
    for (int i = 0; i < max_check_length; i++) {
        char base = std::toupper(r2_seq[i]);
        
        // 距离权重
        double distance_weight = params.calculateDistanceWeight(i);
        
        // 对于R2，GA是目标（CT的complement），CT是惩罚
        bool is_ga = (base == 'G' || base == 'A');
        bool is_ct = (base == 'C' || base == 'T');
        
        if (is_ga) {
            // GA碱基：累加连续计数，应用衰减权重
            consecutive_ga_count++;
            consecutive_ct_count = 0; // 重置CT计数
            
            // GA权重衰减
            current_ga_weight *= params.consecutive_ct_decay_rate;
            current_ct_penalty_multiplier = 1.0; // 重置CT惩罚倍数
            
            // GA碱基加分（按照对应的CT分数）
            double base_score = (base == 'G') ? params.c_score : params.t_score;
            double final_weight = distance_weight * current_ga_weight;
            current_score += base_score * final_weight;
            
            // 濒死模式复活检查：连续3个GA可以复活
            if (in_death_mode && consecutive_ga_count >= 3) {
                in_death_mode = false;
                death_mode_start = -1;
            }
            
        } else if (is_ct) {
            // CT碱基：高惩罚且权重递增
            consecutive_ct_count++;
            consecutive_ga_count = 0; // 重置GA计数
            
            // CT惩罚权重递增
            current_ct_penalty_multiplier *= (1.0 + params.ag_penalty_increment);
            current_ga_weight = 1.0; // 重置GA权重
            
            // CT碱基扣分 - 前8个位置扣分减半
            double base_penalty = params.ag_penalty_base;
            double final_penalty = distance_weight * current_ct_penalty_multiplier;
            
            // 前8个位置扣分减半
            if (i < 10) {
                final_penalty *= 0.5;
            }
            
            current_score += base_penalty * final_penalty;
            
            // 检查是否进入濒死模式
            if (current_score < -20.0 && !in_death_mode) {  // 濒死阈值
                in_death_mode = true;
                death_mode_start = i;
            }
            
        } else {
            // N或其他碱基
            consecutive_ga_count = 0;
            consecutive_ct_count = 0;
            current_ga_weight = 1.0;
            current_ct_penalty_multiplier = 1.0;
            
            double base_score = params.getBaseScore(base);
            current_score += base_score * distance_weight;
        }
        
        // 记录下一个分数
        if (best_cutoff > 0 && !has_next_score) {
            cutoff_next_score = current_score;
            has_next_score = true;
        }
        
        // 无条件更新最佳cutoff（分数没cutoff就一直走）
        int potential_cutoff = i + 1;
        
        // 检查概率和边界条件
        if (potential_cutoff > 0) {
            // 计算要切除的序列概率（R2头部对应CT tail的complement）
            // 将GA序列转换为对应的CT序列计算概率
            std::string cut_sequence = r2_seq.substr(0, potential_cutoff);
            std::string equivalent_ct_sequence = "";
            
            for (char c : cut_sequence) {
                char upper_c = std::toupper(c);
                if (upper_c == 'G') {
                    equivalent_ct_sequence += 'C';
                } else if (upper_c == 'A') {
                    equivalent_ct_sequence += 'T';
                } else {
                    equivalent_ct_sequence += upper_c;
                }
            }
            
            double cut_probability = params.calculateCTTailProbability(equivalent_ct_sequence, 0, equivalent_ct_sequence.length());
            
            // 检查切除边界
            bool valid_boundary = params.isValidCutBoundary(r2_seq, potential_cutoff, false);
            
            // 如果概率满足条件且边界有效，更新cutoff
            if (cut_probability >= params.min_ct_tail_probability && valid_boundary) {
                best_cutoff = potential_cutoff;
                best_cutoff_score = current_score;
            }
        }
        
        // 濒死模式超时检查
        if (in_death_mode && (i - death_mode_start) > 10) {  // 10次后超时
            death_count++;
            if (death_count >= max_death_count) {
                break; // 濒死2次后彻底停止
            }
            // 重置濒死状态，继续下一轮濒死
            in_death_mode = false;
            death_mode_start = -1;
        }
    }
    
    if (best_cutoff > 0) {
        result.r1_trim_pos = 0; // 不trim R1
        result.r2_trim_pos = best_cutoff; // 从R2开头trim
        result.final_score = best_cutoff_score;
        
        // 生成详细信息
        std::string cut_sequence = r2_seq.substr(0, best_cutoff);
        std::string equivalent_ct_sequence = "";
        for (char c : cut_sequence) {
            char upper_c = std::toupper(c);
            if (upper_c == 'G') {
                equivalent_ct_sequence += 'C';
            } else if (upper_c == 'A') {
                equivalent_ct_sequence += 'T';
            } else {
                equivalent_ct_sequence += upper_c;
            }
        }
        double final_probability = params.calculateCTTailProbability(equivalent_ct_sequence, 0, equivalent_ct_sequence.length());
        
        result.score_detail = "R2_score=" + std::to_string(best_cutoff_score) + 
                             ",prob=" + std::to_string(final_probability) +
                             ",len=" + std::to_string(best_cutoff);
        
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
    int check_length_r1 = std::min(cut_length_from_end + 5, (int)original_r1_seq.length());
    if (check_length_r1 > 0) {
        std::string r1_end = original_r1_seq.substr(original_r1_seq.length() - check_length_r1);
        for (char base : r1_end) {
            if (params.isTargetBase(base)) {
                r1_end_c_count++;
            }
        }
    }
    
    // 检查原始R2开头的polyG含量（对应polyC）
    int r2_start_g_count = 0;
    int check_length_r2 = std::min(cut_length_from_end + 5, (int)original_r2_seq.length());
    if (check_length_r2 > 0) {
        std::string r2_start = original_r2_seq.substr(0, check_length_r2);
        for (char base : r2_start) {
            // 对于R2，检查complement bases
            bool is_complement = false;
            if (params.target_bases.find("C") != std::string::npos && std::toupper(base) == 'G') {
                is_complement = true;
            } else if (params.target_bases.find("T") != std::string::npos && std::toupper(base) == 'A') {
                is_complement = true;
            }
            if (is_complement || std::toupper(base) == 'N') {
                r2_start_g_count++;
            }
        }
    }
    
    // 生物学对应关系：
    // 如果R1末尾有polyC且R2开头有polyG，它们在merge后都变成同一段polyC tail
    // 我们需要同时trim它们
    
    double r1_target_ratio = (double)r1_end_c_count / check_length_r1;
    double r2_target_ratio = (double)r2_start_g_count / check_length_r2;
    
    if (r1_target_ratio >= 0.3 && r2_target_ratio >= 0.3) {
        // 两者都有明显的target特征 - 同时trim
        r1_trim = std::min(r1_end_c_count, cut_length_from_end);
        r2_trim = std::min(r2_start_g_count, cut_length_from_end);
    } else if (r1_target_ratio >= 0.3) {
        // 只有R1有polyC - 只trim R1
        r1_trim = std::min(r1_end_c_count + 3, cut_length_from_end); // 稍微aggressive一点
        r2_trim = 0;
    } else if (r2_target_ratio >= 0.3) {
        // 只有R2有polyG - 只trim R2
        r1_trim = 0;
        r2_trim = std::min(r2_start_g_count + 3, cut_length_from_end); // 稍微aggressive一点
    } else {
        // 都没有明显特征 - 但算法检测到了polyC tail，说明确实有问题
        // 按检测到的长度分配，优先考虑更可能的来源
        if (r1_end_c_count >= r2_start_g_count) {
            r1_trim = std::min(cut_length_from_end, check_length_r1);
            r2_trim = std::max(0, cut_length_from_end - r1_trim);
        } else {
            r2_trim = std::min(cut_length_from_end, check_length_r2);
            r1_trim = std::max(0, cut_length_from_end - r2_trim);
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
    
    std::vector<int> window_sizes = {25, 30, 40, 50};
    std::vector<double> initial_scores = {-10.0, -15.0, -20.0};
    
    // PolyC配置矩阵 - 更aggressive参数
    std::vector<double> c_scores = {15.0, 18.0, 20.0, 25.0};
    std::vector<double> penalty_scores = {-30.0, -35.0, -40.0, -45.0};
    std::vector<double> cutoff_thresholds = {-5.0, -8.0, -10.0, -12.0};
    
    for (int ws : window_sizes) {
        for (double c_score : c_scores) {
            for (double penalty : penalty_scores) {
                for (double init_score : initial_scores) {
                    for (double cutoff_thresh : cutoff_thresholds) {
                        TrimParams params = TrimParams::createPolyCConfig(ws);
                        params.c_score = c_score;
                        params.a_score = penalty;
                        params.g_score = penalty;
                        params.t_score = penalty;
                        params.initial_score = init_score;
                        params.cutoff_threshold = cutoff_thresh;
                        param_list.push_back(params);
                    }
                }
            }
        }
    }
    
    // CT tail配置矩阵
    std::vector<double> t_scores = {12.0, 15.0, 18.0}; // T分数通常比C分数低
    std::vector<double> ct_ratios = {0.7, 0.8, 0.9, 1.0}; // T/C 比例
    
    for (int ws : window_sizes) {
        for (double c_score : c_scores) {
            for (double t_ratio : ct_ratios) {
                for (double penalty : penalty_scores) {
                    for (double init_score : initial_scores) {
                        for (double cutoff_thresh : cutoff_thresholds) {
                            TrimParams params = TrimParams::createCTConfig(ws);
                            params.c_score = c_score;
                            params.t_score = c_score * t_ratio;
                            params.a_score = penalty;
                            params.g_score = penalty;
                            params.initial_score = init_score;
                            params.cutoff_threshold = cutoff_thresh;
                            param_list.push_back(params);
                        }
                    }
                }
            }
        }
    }
    
    return param_list;
}