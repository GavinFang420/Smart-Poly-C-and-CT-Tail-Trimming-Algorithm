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
    int window_size = 30;
    
    // 碱基权重配置 - 支持ACGT四种碱基独立配置
    double a_score = -20.0;        // A的分数（通常是惩罚）
    double c_score = 10.0;         // C的分数（polyC tail target）
    double g_score = -20.0;        // G的分数（通常是惩罚，除非检测polyG）
    double t_score = -20.0;        // T的分数（通常是惩罚）
    double n_score = 5.0;          // N的分数（通常当作target处理）
    
    // 衰减函数配置
    enum DecayFunction {
        LINEAR,      // 线性衰减: weight = 1.0 - decay_rate * position / window_size
        EXPONENTIAL, // 指数衰减: weight = exp(-decay_rate * position / window_size)
        QUADRATIC,   // 二次衰减: weight = 1.0 - decay_rate * (position/window_size)^2
        LOGARITHMIC  // 对数衰减: weight = 1.0 - decay_rate * log(1 + position) / log(1 + window_size)
    };
    
    DecayFunction distance_decay_func = LINEAR;
    double distance_decay_rate = 0.95;     // 距离衰减率
    double min_distance_weight = 0.05;     // 最小距离权重
    
    DecayFunction consecutive_decay_func = EXPONENTIAL;
    double consecutive_decay_rate = 2.0;   // 连续衰减率 (对于指数: 1/2^n)
    double min_consecutive_weight = 0.01;  // 最小连续权重
    
    double initial_score = 0.0;           // 初始分数
    
    // 目标检测配置
    std::string target_bases = "C";        // 目标碱基序列 "C", "CT", "CG" 等
    bool treat_n_as_target = true;         // 是否将N当作目标碱基处理
    
    // 构造函数
    TrimParams(int ws = 30, double init_score = 0.0) 
        : window_size(ws), initial_score(init_score) {
        // 默认配置为polyC检测
    }
    
    // 预设配置函数
    static TrimParams createPolyCConfig(int window_size = 30) {
        TrimParams params(window_size, 0.0);
        params.c_score = 10.0;
        params.target_bases = "C";
        return params;
    }
    
    static TrimParams createCTConfig(int window_size = 30) {
        TrimParams params(window_size, 0.0);
        params.c_score = 10.0;
        params.t_score = 8.0;          // CT tail中T的分数略低于C
        params.target_bases = "CT";
        params.consecutive_decay_rate = 1.5;  // CT tail连续衰减更温和
        return params;
    }
    
    static TrimParams createPolyGConfig(int window_size = 30) {
        TrimParams params(window_size, 0.0);
        params.g_score = 10.0;
        params.target_bases = "G";
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
};

struct TrimResult {
    int r1_trim_pos;    // Position to trim R1 (from end)
    int r2_trim_pos;    // Position to trim R2 (from start)
    double final_score;
    std::string score_detail; // 详细的score信息 "正数 -> 负数"
    bool is_valid;
    
    TrimResult() : r1_trim_pos(0), r2_trim_pos(0), final_score(0.0), is_valid(false) {}
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
};

#endif
