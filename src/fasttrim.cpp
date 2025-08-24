#include "fasttrim.h"
#include <algorithm>
#include <iostream>

// R1尾部CT tail检测 (无濒死，简单模式)
int FastTrimmer::analyzeR1CTTail(const std::string& r1_seq) {
    if (r1_seq.empty()) {
        return 0;
    }
    
    double current_score = params.initial_score;
    int best_cutoff = 0;
    
    int consecutive_ct_count = 0;
    int consecutive_ag_count = 0;
    double ct_weight = 1.0;
    double ag_penalty_multiplier = 1.0;
    
    // 从尾部开始检测，最多检测20bp
    int max_check = std::min(params.r1_max_trim_length, (int)r1_seq.length());
    
    for (int i = 0; i < max_check; i++) {
        int pos = r1_seq.length() - 1 - i; // 从末尾开始
        char base = std::toupper(r1_seq[pos]);
        
        bool is_ct = (base == 'C' || base == 'T');
        bool is_ag = (base == 'A' || base == 'G');
        
        if (is_ct || base == 'N') {
            // CT或N碱基加分
            consecutive_ct_count++;
            consecutive_ag_count = 0;
            
            // 连续加分衰减
            ct_weight *= params.r1_consecutive_decay;
            ag_penalty_multiplier = 1.0; // 重置扣分倍率
            
            double base_score = (base == 'N') ? params.r1_n_score : 
                               (base == 'C') ? params.r1_c_score : params.r1_t_score;
            current_score += base_score * ct_weight;
            
        } else if (is_ag) {
            // AG碱基扣分
            consecutive_ag_count++;
            consecutive_ct_count = 0;
            
            // 连续扣分增加倍率
            ag_penalty_multiplier *= params.r1_penalty_increase;
            ct_weight = 1.0; // 重置加分权重
            
            current_score += params.r1_penalty * ag_penalty_multiplier;
        } else {
            // 其他碱基
            consecutive_ct_count = 0;
            consecutive_ag_count = 0;
            ct_weight = 1.0;
            ag_penalty_multiplier = 1.0;
            
            current_score += params.r1_penalty; // 默认扣分
        }
        
        // 只要分数为正就更新cutoff (无濒死限制)
        if (current_score > 0.0) {
            best_cutoff = i + 1;
        }
    }
    
    return best_cutoff;
}

// R2头部GA tail检测 (有濒死复活)
int FastTrimmer::analyzeR2GATail(const std::string& r2_seq) {
    if (r2_seq.empty()) {
        return 0;
    }
    
    double current_score = params.initial_score;
    int best_cutoff = 0;
    
    // 濒死状态控制
    bool in_death_mode = false;
    int death_mode_start = -1;
    int death_count = 0;
    
    int consecutive_ga_count = 0;
    int consecutive_ct_count = 0;
    double ga_weight = 1.0;
    double ct_penalty_multiplier = 1.0;
    
    // 从头部开始检测，最多检测20bp
    int max_check = std::min(params.r2_max_trim_length, (int)r2_seq.length());
    
    for (int i = 0; i < max_check; i++) {
        char base = std::toupper(r2_seq[i]);
        
        bool is_ga = (base == 'G' || base == 'A');
        bool is_ct = (base == 'C' || base == 'T');
        
        // 前6个碱基免伤80%
        double immunity_factor = (i < 6) ? params.r2_immunity_ratio : 1.0;
        
        if (is_ga || base == 'N') {
            // GA或N碱基加分
            consecutive_ga_count++;
            consecutive_ct_count = 0;
            
            // 连续加分衰减
            ga_weight *= params.r2_consecutive_decay;
            ct_penalty_multiplier = 1.0; // 重置扣分倍率
            
            double base_score = (base == 'N') ? params.r2_n_score : 
                               (base == 'G') ? params.r2_g_score : params.r2_a_score;
            current_score += base_score * ga_weight;
            
            // 濒死复活检查：连续3个GA可以复活
            if (in_death_mode && consecutive_ga_count >= 3) {
                in_death_mode = false;
                death_mode_start = -1;
                current_score = 0.0; // 复活时重置分数到0
            }
            
        } else if (is_ct) {
            // CT碱基扣分
            consecutive_ct_count++;
            consecutive_ga_count = 0;
            
            // 连续扣分增加倍率
            ct_penalty_multiplier *= params.r2_penalty_increase;
            ga_weight = 1.0; // 重置加分权重
            
            // 应用免伤和扣分
            double penalty = params.r2_penalty * ct_penalty_multiplier * immunity_factor;
            current_score += penalty;
            
            // 检查是否进入濒死模式
            if (current_score < -10.0 && !in_death_mode) {
                in_death_mode = true;
                death_mode_start = i;
            }
            
        } else {
            // 其他碱基
            consecutive_ga_count = 0;
            consecutive_ct_count = 0;
            ga_weight = 1.0;
            ct_penalty_multiplier = 1.0;
            
            // 应用免伤和默认扣分
            double penalty = params.r2_penalty * immunity_factor;
            current_score += penalty;
            
            if (current_score < -10.0 && !in_death_mode) {
                in_death_mode = true;
                death_mode_start = i;
            }
        }
        
        // 更新最佳cutoff
        if (current_score > 0.0) {
            best_cutoff = i + 1;
        }
        
        // 濒死模式超时检查
        if (in_death_mode && (i - death_mode_start) >= params.r2_death_timeout) {
            death_count++;
            if (death_count >= params.r2_death_count_limit) {
                break; // 濒死2次后停止
            }
            // 重置濒死状态，重置分数到0
            in_death_mode = false;
            death_mode_start = -1;
            current_score = 0.0;
        }
    }
    
    return best_cutoff;
}

// 主要处理函数
FastTrimResult FastTrimmer::trimReadPair(const std::string& r1_seq, const std::string& r2_seq) {
    FastTrimResult result;
    
    // 分别分析R1和R2
    result.r1_trim_pos = analyzeR1CTTail(r1_seq);
    result.r2_trim_pos = analyzeR2GATail(r2_seq);
    
    // 设置有效性
    result.r1_valid = (result.r1_trim_pos > 0);
    result.r2_valid = (result.r2_trim_pos > 0);
    
    // 记录分数（简化版，主要记录trim长度）
    result.r1_score = result.r1_trim_pos;
    result.r2_score = result.r2_trim_pos;
    
    return result;
}

// 应用trim结果
std::pair<std::string, std::string> FastTrimmer::applyTrim(const std::string& r1_seq, 
                                                          const std::string& r2_seq, 
                                                          const FastTrimResult& result) {
    std::string trimmed_r1 = r1_seq;
    std::string trimmed_r2 = r2_seq;
    
    // 从R1尾部trim
    if (result.r1_valid && result.r1_trim_pos > 0 && result.r1_trim_pos < (int)r1_seq.length()) {
        trimmed_r1 = r1_seq.substr(0, r1_seq.length() - result.r1_trim_pos);
    }
    
    // 从R2头部trim
    if (result.r2_valid && result.r2_trim_pos > 0 && result.r2_trim_pos < (int)r2_seq.length()) {
        trimmed_r2 = r2_seq.substr(result.r2_trim_pos);
    }
    
    return std::make_pair(trimmed_r1, trimmed_r2);
}

// 预设配置
FastTrimParams FastTrimmer::createCTConfig() {
    FastTrimParams params;
    
    // R1 CT tail检测参数
    params.r1_c_score = 10.0;
    params.r1_t_score = 10.0;
    params.r1_n_score = 10.0;
    params.r1_penalty = -20.0;
    params.r1_max_trim_length = 20;
    params.r1_consecutive_decay = 0.9;
    params.r1_penalty_increase = 1.2;
    
    // R2 GA tail检测参数  
    params.r2_g_score = 10.0;
    params.r2_a_score = 10.0;
    params.r2_n_score = 10.0;
    params.r2_penalty = -20.0;
    params.r2_max_trim_length = 20;
    params.r2_consecutive_decay = 0.9;
    params.r2_penalty_increase = 1.2;
    params.r2_death_count_limit = 2;
    params.r2_death_timeout = 10;
    params.r2_immunity_ratio = 0.2; // 前6个碱基80%免伤 = 20%伤害
    
    params.initial_score = 0.0;
    
    return params;
}

FastTrimParams FastTrimmer::createPolyCConfig() {
    FastTrimParams params;
    
    // R1 C tail检测参数 (T变成扣分)
    params.r1_c_score = 10.0;
    params.r1_t_score = -20.0;   // 纯C模式下T扣分
    params.r1_n_score = 10.0;
    params.r1_penalty = -20.0;
    params.r1_max_trim_length = 20;
    params.r1_consecutive_decay = 0.9;
    params.r1_penalty_increase = 1.2;
    
    // R2 G tail检测参数 (A变成扣分)
    params.r2_g_score = 10.0;
    params.r2_a_score = -20.0;   // 纯G模式下A扣分
    params.r2_n_score = 10.0;
    params.r2_penalty = -20.0;
    params.r2_max_trim_length = 20;
    params.r2_consecutive_decay = 0.9;
    params.r2_penalty_increase = 1.2;
    params.r2_death_count_limit = 2;
    params.r2_death_timeout = 10;
    params.r2_immunity_ratio = 0.2;
    
    params.initial_score = 0.0;
    
    return params;
}