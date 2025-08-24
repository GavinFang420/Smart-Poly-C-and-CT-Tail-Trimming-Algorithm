/*
MIT License

Copyright (c) 2017 OpenGene
Copyright (c) 2025 GavinFang420 - Smart Poly-C/CT Tail Trimming Algorithm

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SMARTTRIM_H
#define SMARTTRIM_H

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

// Forward declaration
struct MergeResult;

struct TrimParams {
    int window_size = 40;  // 增加默认窗口大小
    
    // 碱基权重配置 - 更aggressive的参数
    double a_score = -35.0;        // A的分数（增加惩罚）
    double c_score = 18.0;         // C的分数（增加奖励）
    double g_score = -35.0;        // G的分数（增加惩罚）
    double t_score = -35.0;        // T的分数（增加惩罚）
    double n_score = 8.0;          // N的分数（轻微增加）
    
    // 衰减函数配置
    enum DecayFunction {
        LINEAR,      // 线性衰减: weight = 1.0 - decay_rate * position / window_size
        EXPONENTIAL, // 指数衰减: weight = exp(-decay_rate * position / window_size)
        QUADRATIC,   // 二次衰减: weight = 1.0 - decay_rate * (position/window_size)^2
        LOGARITHMIC  // 对数衰减: weight = 1.0 - decay_rate * log(1 + position) / log(1 + window_size)
    };
    
    DecayFunction distance_decay_func = LINEAR;
    double distance_decay_rate = 0.75;     // 降低距离衰减率，让远端也能有效trim
    double min_distance_weight = 0.15;     // 提高最小距离权重
    
    DecayFunction consecutive_decay_func = EXPONENTIAL;
    double consecutive_decay_rate = 1.2;   // 降低连续衰减率
    double min_consecutive_weight = 0.05;  // 提高最小连续权重
    
    double initial_score = -15.0;          // 降低初始分数，更容易触发trim
    
    // 目标检测配置
    std::string target_bases = "C";        // 目标碱基序列
    bool treat_n_as_target = true;         // 是否将N当作目标碱基处理
    
    // 新增：概率计算参数
    double ct_c_probability = 0.8;        // CT tail中C的概率
    double ct_t_probability = 0.2;        // CT tail中T的概率  
    double ct_other_probability = 0.1;    // CT tail中其他碱基(非N)的概率
    double ct_n_probability = 1.0;        // N按照1.0概率处理
    double min_ct_tail_probability = 0.005; // 最低CT tail概率阈值(0.5%)
    
    // 权重衰减和惩罚参数
    double consecutive_ct_decay_rate = 0.9;   // 连续CT权重衰减率
    double ag_penalty_base = -50.0;           // AG基础惩罚
    double ag_penalty_increment = 0.2;        // AG惩罚递增率(每次+20%)
    
    // 切除边界检查
    bool require_non_ct_after_cut = true;    // 要求切除位置后必须不是C/T
    
    // 保留的旧参数（为了兼容性）
    double cutoff_threshold = -8.0;        // cutoff阈值（负数表示允许负分数cutoff）
    double force_trim_threshold = 0.45;    // 强制trim阈值（target base比例）
    int max_force_trim_length = 20;        // 最大强制trim长度
    double penalty_multiplier_base = 2.0;  // 惩罚倍数基础值
    double penalty_multiplier_max = 4.0;   // 惩罚倍数最大值
    int death_mode_recovery_bases = 2;     // death mode复活需要的连续target bases
    int death_mode_timeout = 10;           // death mode超时长度
    
    // 构造函数
    TrimParams(int ws = 40, double init_score = -15.0) 
        : window_size(ws), initial_score(init_score) {
        // 默认配置为aggressive polyC检测
    }
    
    // 预设配置函数 - 更aggressive的参数
    static TrimParams createPolyCConfig(int window_size = 40) {
        TrimParams params(window_size, -15.0);
        params.c_score = 20.0;               // 大幅提高C的奖励
        params.a_score = -40.0;              // 增加非C碱基的惩罚
        params.g_score = -40.0;
        params.t_score = -40.0;
        params.target_bases = "C";
        params.distance_decay_rate = 0.7;    
        params.min_distance_weight = 0.2;    
        // C tail概率参数（纯C更严格）
        params.ct_c_probability = 0.95;      // 纯C tail中C概率更高
        params.ct_t_probability = 0.0;       // 纯C tail中不应该有T
        params.ct_other_probability = 0.05;  // 其他碱基概率更低
        params.min_ct_tail_probability = 0.01; // 更高的概率阈值(1%)
        // 权重和惩罚参数
        params.consecutive_ct_decay_rate = 0.95; // 纯C衰减更慢
        params.ag_penalty_base = -60.0;      // 更重的AG惩罚
        params.ag_penalty_increment = 0.25;  // 更快的惩罚递增
        params.require_non_ct_after_cut = false; // 纯C可以更宽松
        return params;
    }
    
    static TrimParams createCTConfig(int window_size = 40) {
        TrimParams params(window_size, -12.0);
        params.c_score = 18.0;               // 高C奖励
        params.t_score = 15.0;               // 高T奖励（略低于C）
        params.a_score = -35.0;              // 增加A/G惩罚
        params.g_score = -35.0;
        params.target_bases = "CT";
        params.distance_decay_rate = 0.75;   
        params.min_distance_weight = 0.18;   
        // CT概率参数
        params.ct_c_probability = 0.8;
        params.ct_t_probability = 0.2;
        params.ct_other_probability = 0.1;
        params.min_ct_tail_probability = 0.005;
        // 权重和惩罚参数
        params.consecutive_ct_decay_rate = 0.9;
        params.ag_penalty_base = -50.0;
        params.ag_penalty_increment = 0.2;
        params.require_non_ct_after_cut = true;
        return params;
    }
    
    static TrimParams createPolyGConfig(int window_size = 40) {
        TrimParams params(window_size, -12.0);
        params.g_score = 18.0;               // 高G奖励
        params.a_score = -35.0;              // 增加非G碱基惩罚
        params.c_score = -35.0;
        params.t_score = -35.0;
        params.target_bases = "G";
        params.consecutive_decay_rate = 1.1;
        params.distance_decay_rate = 0.75;
        params.min_distance_weight = 0.18;
        params.cutoff_threshold = -8.0;
        params.force_trim_threshold = 0.4;
        params.penalty_multiplier_base = 2.2;
        params.penalty_multiplier_max = 4.5;
        return params;
    }
    
    // 获取指定碱基的分数
    double getBaseScore(char base) const {
        switch (std::toupper(base)) {
            case 'A': return a_score;
            case 'C': return c_score;
            case 'G': return g_score;
            case 'T': return t_score;
            case 'N': return n_score;
            default: return a_score; // 默认按A处理
        }
    }
    
    // 检查碱基是否是目标碱基
    bool isTargetBase(char base) const {
        if (treat_n_as_target && std::toupper(base) == 'N') {
            return true;
        }
        return target_bases.find(std::toupper(base)) != std::string::npos;
    }
    
    // 计算距离权重
    double calculateDistanceWeight(int position) const {
        double ratio = (double)position / window_size;
        double weight = 1.0;
        
        switch (distance_decay_func) {
            case LINEAR:
                weight = 1.0 - distance_decay_rate * ratio;
                break;
            case EXPONENTIAL:
                weight = std::exp(-distance_decay_rate * ratio);
                break;
            case QUADRATIC:
                weight = 1.0 - distance_decay_rate * ratio * ratio;
                break;
            case LOGARITHMIC:
                weight = 1.0 - distance_decay_rate * std::log(1.0 + position) / std::log(1.0 + window_size);
                break;
        }
        
        return std::max(min_distance_weight, weight);
    }
    
    // 计算连续权重
    double calculateConsecutiveWeight(int consecutive_count) const {
        if (consecutive_count <= 1) return 1.0;
        
        double weight = 1.0;
        switch (consecutive_decay_func) {
            case LINEAR:
                weight = 1.0 - distance_decay_rate * (consecutive_count - 1) / 10.0;
                break;
            case EXPONENTIAL:
                weight = 1.0 / std::pow(consecutive_decay_rate, consecutive_count - 1);
                break;
            case QUADRATIC:
                weight = 1.0 / std::pow(consecutive_count, consecutive_decay_rate);
                break;
            case LOGARITHMIC:
                weight = 1.0 / std::log(consecutive_decay_rate + consecutive_count - 1);
                break;
        }
        
        return std::max(min_consecutive_weight, weight);
    }
    
    // 计算惩罚倍数
    double calculatePenaltyMultiplier(double distance_weight) const {
        return penalty_multiplier_base + (penalty_multiplier_max - penalty_multiplier_base) * distance_weight;
    }
    
    // 新增：CT tail概率计算函数
    double calculateCTTailProbability(const std::string& sequence, int start_pos, int end_pos) const {
        if (start_pos >= end_pos || start_pos < 0 || end_pos > (int)sequence.length()) {
            return 0.0;
        }
        
        double total_probability = 1.0;
        int length = end_pos - start_pos;
        
        for (int i = start_pos; i < end_pos; i++) {
            char base = std::toupper(sequence[i]);
            double base_prob = 0.0;
            
            switch (base) {
                case 'C':
                    base_prob = ct_c_probability;
                    break;
                case 'T':
                    base_prob = ct_t_probability;
                    break;
                case 'N':
                    base_prob = ct_n_probability;
                    break;
                default:
                    base_prob = ct_other_probability;
                    break;
            }
            
            total_probability *= base_prob;
        }
        
        return total_probability;
    }
    
    // 检查切除边界是否有效
    bool isValidCutBoundary(const std::string& sequence, int cut_pos, bool is_tail_cut = true) const {
        if (!require_non_ct_after_cut) {
            return true;
        }
        
        if (is_tail_cut) {
            // 尾部切除：检查cut_pos之前的位置（即保留的最后一个碱基之后）
            if (cut_pos >= (int)sequence.length()) {
                return true; // 切到头了，没问题
            }
            int check_pos = sequence.length() - cut_pos;
            if (check_pos < (int)sequence.length()) {
                char base = std::toupper(sequence[check_pos]);
                return (base != 'C' && base != 'T');
            }
        } else {
            // 头部切除：检查cut_pos位置（即保留的第一个碱基）
            if (cut_pos >= (int)sequence.length()) {
                return true;
            }
            char base = std::toupper(sequence[cut_pos]);
            return (base != 'C' && base != 'T');
        }
        
        return true;
    }
};

struct TrimResult {
    int r1_trim_pos;    // Position to trim R1 (from end)
    int r2_trim_pos;    // Position to trim R2 (from start)
    double final_score;
    std::string score_detail; // 详细的score信息 "正数 -> 负数"
    bool is_valid;
    bool force_trimmed; // 新增：标记是否为强制trim
    
    TrimResult() : r1_trim_pos(0), r2_trim_pos(0), final_score(0.0), 
                   is_valid(false), force_trimmed(false) {}
};

class SmartTrimmer {
private:
    TrimParams params;
    
    // 新的动态递进检测算法 - 核心函数
    TrimResult analyzeMergedSequence(const MergeResult& merge_result,
                                   const std::string& original_r1_seq,
                                   const std::string& original_r2_seq);
    
    TrimResult analyzeR2PolyGTail(const std::string& r2_seq);
    
    std::pair<int, int> mapMergedPositionToOriginal(int cut_length_from_end, 
                                                   const MergeResult& merge_result,
                                                   const std::string& original_r1_seq,
                                                   const std::string& original_r2_seq);
    
    // 新增：预检查函数
    double performPreCheck(const std::string& sequence, bool is_tail_check = true) const;
    
    // 新增：强制trim检查
    int checkForceTrimmingCondition(const std::string& sequence, 
                                   double target_ratio, bool is_tail_check = true) const;

public:
    SmartTrimmer(const TrimParams& p) : params(p) {}
    
    // Main function to find optimal trim positions
    // 新算法: 动态递进检测，支持濒死复活机制，支持ACGT独立权重配置
    TrimResult findOptimalTrimPositions(
        const std::string& r1_seq, 
        const std::string& r2_seq
    );
    
    // Trim the reads based on result
    std::pair<std::string, std::string> trimReads(
        const std::string& r1_seq,
        const std::string& r2_seq,
        const TrimResult& result
    );
    
    // Generate multiple parameter configurations for testing
    static std::vector<TrimParams> generateParameterMatrix();
    
    // 新增：获取参数配置
    const TrimParams& getParams() const { return params; }
    void setParams(const TrimParams& p) { params = p; }
};

#endif