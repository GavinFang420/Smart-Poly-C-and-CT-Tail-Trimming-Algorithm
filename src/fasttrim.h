#ifndef FASTTRIM_H
#define FASTTRIM_H

#include <string>
#include <vector>

struct FastTrimParams {
    // R1 参数 (CT tail检测)
    double r1_c_score = 10.0;           // C加分
    double r1_t_score = 10.0;           // T加分  
    double r1_n_score = 10.0;           // N加分
    double r1_penalty = -20.0;          // A/G扣分
    int r1_max_trim_length = 20;        // 最大trim长度
    double r1_consecutive_decay = 0.9;  // 连续加分衰减
    double r1_penalty_increase = 1.2;   // 连续扣分增加倍率
    
    // R2 参数 (GA tail检测，有濒死复活)
    double r2_g_score = 10.0;           // G加分
    double r2_a_score = 10.0;           // A加分
    double r2_n_score = 10.0;           // N加分  
    double r2_penalty = -20.0;          // C/T扣分
    int r2_max_trim_length = 20;        // 最大trim长度
    double r2_consecutive_decay = 0.9;  // 连续加分衰减
    double r2_penalty_increase = 1.2;   // 连续扣分增加倍率
    int r2_death_count_limit = 2;       // 濒死次数限制
    int r2_death_timeout = 10;          // 濒死超时长度
    double r2_immunity_ratio = 0.8;     // 前6个碱基免伤比例
    
    // 通用参数
    double initial_score = 0.0;         // 初始分数
};

struct FastTrimResult {
    int r1_trim_pos;     // R1从尾部trim的长度
    int r2_trim_pos;     // R2从头部trim的长度
    double r1_score;     // R1最终分数
    double r2_score;     // R2最终分数
    bool r1_valid;       // R1是否有效trim
    bool r2_valid;       // R2是否有效trim
    
    FastTrimResult() : r1_trim_pos(0), r2_trim_pos(0), 
                      r1_score(0.0), r2_score(0.0), 
                      r1_valid(false), r2_valid(false) {}
};

class FastTrimmer {
private:
    FastTrimParams params;
    
    // R1尾部CT tail检测 (无濒死，简单模式)
    int analyzeR1CTTail(const std::string& r1_seq);
    
    // R2头部GA tail检测 (有濒死复活)
    int analyzeR2GATail(const std::string& r2_seq);

public:
    FastTrimmer(const FastTrimParams& p = FastTrimParams()) : params(p) {}
    
    // 主要处理函数
    FastTrimResult trimReadPair(const std::string& r1_seq, const std::string& r2_seq);
    
    // 应用trim结果
    std::pair<std::string, std::string> applyTrim(const std::string& r1_seq, 
                                                  const std::string& r2_seq, 
                                                  const FastTrimResult& result);
    
    // 参数设置
    void setParams(const FastTrimParams& p) { params = p; }
    const FastTrimParams& getParams() const { return params; }
    
    // 预设配置
    static FastTrimParams createCTConfig();
    static FastTrimParams createPolyCConfig();
};

#endif