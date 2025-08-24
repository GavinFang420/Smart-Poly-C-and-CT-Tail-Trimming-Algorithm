#ifndef CONSTANTTRIM_H
#define CONSTANTTRIM_H

#include <string>
#include <vector>
#include <map>

struct ConstantTrimParams {
    // 基本trim参数
    int r1_trim_length = 15;        // R1从尾部切除的固定长度
    int r2_trim_length = 15;        // R2从头部切除的固定长度
    int min_remaining_length = 50;   // 最少保留的序列长度
    
    // 质量控制参数
    bool use_quality_check = false;  // 是否使用质量检查来调整trim长度
    double high_quality_threshold = 30.0;  // 高质量阈值（Phred score）
    double low_quality_threshold = 10.0;   // 低质量阈值（Phred score）
    int quality_adjustment = 3;      // 基于质量的trim长度调整量
};

struct BaseComposition {
    int a_count = 0;
    int t_count = 0;
    int c_count = 0;
    int g_count = 0;
    int n_count = 0;
    int other_count = 0;
    int total_bases = 0;
    
    double a_percent = 0.0;
    double t_percent = 0.0;
    double c_percent = 0.0;
    double g_percent = 0.0;
    double n_percent = 0.0;
    double gc_percent = 0.0;
};

struct ConstantTrimResult {
    int r1_trim_pos;                // R1从尾部trim的长度
    int r2_trim_pos;                // R2从头部trim的长度
    bool r1_valid;                  // R1是否有效trim
    bool r2_valid;                  // R2是否有效trim
    
    // 剩余长度统计
    int r1_remaining_length = 0;    // R1剩余长度
    int r2_remaining_length = 0;    // R2剩余长度
    
    // Trim效率统计
    double r1_trim_efficiency = 0.0; // R1 trim效率百分比
    double r2_trim_efficiency = 0.0; // R2 trim效率百分比
    
    // Trim区域碱基组成分析
    BaseComposition r1_trim_composition;
    BaseComposition r2_trim_composition;
    
    // 用于批量测试时记录测试的trim长度
    int tested_trim_length = 0;
    
    ConstantTrimResult() : r1_trim_pos(0), r2_trim_pos(0), 
                          r1_valid(false), r2_valid(false) {}
};

// 用于分析不同trim长度效果的统计结构
struct TrimLengthStats {
    int trim_length = 0;
    int total_reads = 0;
    
    // R1统计
    int r1_trimmed_count = 0;
    double r1_trim_rate = 0.0;          // trim的read百分比
    double avg_r1_trim_length = 0.0;     // 平均trim长度
    double avg_remaining_r1_length = 0.0; // 平均剩余长度
    double avg_r1_gc_in_trimmed = 0.0;   // trim区域平均GC含量
    
    // R2统计
    int r2_trimmed_count = 0;
    double r2_trim_rate = 0.0;          // trim的read百分比
    double avg_r2_trim_length = 0.0;     // 平均trim长度
    double avg_remaining_r2_length = 0.0; // 平均剩余长度
    double avg_r2_gc_in_trimmed = 0.0;   // trim区域平均GC含量
    
    // 内部累计统计（用于计算平均值）
    int total_r1_bases_trimmed = 0;
    int total_r2_bases_trimmed = 0;
    int total_remaining_r1_bases = 0;
    int total_remaining_r2_bases = 0;
    double total_r1_gc_in_trimmed = 0.0;
    double total_r2_gc_in_trimmed = 0.0;
};

struct TrimLengthAnalysis {
    std::map<int, TrimLengthStats> trim_stats; // key: trim_length, value: stats
};

class ConstantTrimmer {
private:
    ConstantTrimParams params;
    
    // R1尾部固定切除
    int analyzeR1ConstantTrim(const std::string& r1_seq);
    
    // R2头部固定切除
    int analyzeR2ConstantTrim(const std::string& r2_seq);
    
    // 基于质量的固定切除
    int analyzeR1ConstantTrimWithQuality(const std::string& r1_seq, const std::string& r1_qual);
    int analyzeR2ConstantTrimWithQuality(const std::string& r2_seq, const std::string& r2_qual);
    
    // 分析trim区域碱基组成
    BaseComposition analyzeTrimRegionComposition(const std::string& seq, int start, int end);

public:
    ConstantTrimmer(const ConstantTrimParams& p = ConstantTrimParams()) : params(p) {}
    
    // 主要处理函数 - 仅序列
    ConstantTrimResult trimReadPair(const std::string& r1_seq, const std::string& r2_seq);
    
    // 主要处理函数 - 序列和质量
    ConstantTrimResult trimReadPair(const std::string& r1_seq, const std::string& r2_seq,
                                   const std::string& r1_qual, const std::string& r2_qual);
    
    // 应用trim结果
    std::pair<std::string, std::string> applyTrim(const std::string& r1_seq, 
                                                  const std::string& r2_seq, 
                                                  const ConstantTrimResult& result);
    
    // 应用trim结果到质量字符串
    std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string>> 
    applyTrimWithQuality(const std::string& r1_seq, const std::string& r2_seq,
                        const std::string& r1_qual, const std::string& r2_qual,
                        const ConstantTrimResult& result);
    
    // 参数设置
    void setParams(const ConstantTrimParams& p) { params = p; }
    const ConstantTrimParams& getParams() const { return params; }
    
    // 批量测试不同的trim长度 - 针对单个read pair
    std::vector<ConstantTrimResult> testTrimLengthRange(const std::string& r1_seq, 
                                                       const std::string& r2_seq,
                                                       const std::string& r1_qual = "",
                                                       const std::string& r2_qual = "",
                                                       int min_trim = 10, int max_trim = 18);
    
    // 分析不同trim长度对序列质量的影响 - 针对多个read pairs
    TrimLengthAnalysis analyzeTrimLengthImpact(const std::vector<std::string>& r1_seqs,
                                              const std::vector<std::string>& r2_seqs,
                                              const std::vector<std::string>& r1_quals = {},
                                              const std::vector<std::string>& r2_quals = {},
                                              int min_trim = 10, int max_trim = 18);
    
    // 预设配置
    static ConstantTrimParams createDefaultConfig();
    static ConstantTrimParams createR1TailConfig(int trim_length);
    static ConstantTrimParams createR2HeadConfig(int trim_length);
    static ConstantTrimParams createQualityAwareConfig(int trim_length);
    
    // 工具函数
    static bool validateSequence(const std::string& seq);
    static double calculateSequenceGC(const std::string& seq);
    static BaseComposition analyzeSequenceComposition(const std::string& seq);
    
    // 比较不同trim长度的效果
    static std::vector<int> recommendOptimalTrimLengths(const TrimLengthAnalysis& analysis,
                                                       double min_remaining_length_ratio = 0.7);
};

// 工具命名空间
namespace ConstantTrimUtils {
    // 序列验证
    bool isValidDNASequence(const std::string& seq);
    bool isValidQualityString(const std::string& qual);
    
    // 质量分析
    double calculateAverageQuality(const std::string& qual);
    double calculateRegionAverageQuality(const std::string& qual, int start, int end);
    
    // 碱基组成分析
    BaseComposition getRegionComposition(const std::string& seq, int start, int end);
    double calculateGCContent(const BaseComposition& comp);
    
    // 统计分析
    double calculateTrimEffectiveness(const TrimLengthStats& stats);
    int findOptimalTrimLength(const std::map<int, TrimLengthStats>& stats_map, 
                             const std::string& criterion = "balance"); // "balance", "aggressive", "conservative"
    
    // 可视化数据准备
    std::vector<double> extractMetricSeries(const std::map<int, TrimLengthStats>& stats_map,
                                           const std::string& metric); // "r1_trim_rate", "r2_trim_rate", "avg_remaining_length", etc.
    
    // 序列比较
    double calculateSequenceSimilarity(const std::string& seq1, const std::string& seq2);
    int calculateEditDistance(const std::string& seq1, const std::string& seq2);
    
    // R1尾部和R2头部特异性分析
    struct R1TailAnalysis {
        double ct_ratio = 0.0;          // C+T占比
        double poly_c_length = 0.0;     // 平均polyC长度
        double poly_t_length = 0.0;     // 平均polyT长度
        double n_ratio = 0.0;           // N碱基占比
        double avg_quality = 0.0;       // 平均质量
    };
    
    struct R2HeadAnalysis {
        double ga_ratio = 0.0;          // G+A占比
        double poly_g_length = 0.0;     // 平均polyG长度
        double poly_a_length = 0.0;     // 平均polyA长度
        double n_ratio = 0.0;           // N碱基占比
        double avg_quality = 0.0;       // 平均质量
    };
    
    R1TailAnalysis analyzeR1TailRegion(const std::string& seq, const std::string& qual, 
                                      int tail_length);
    R2HeadAnalysis analyzeR2HeadRegion(const std::string& seq, const std::string& qual, 
                                      int head_length);
    
    // 批量分析工具
    std::vector<R1TailAnalysis> batchAnalyzeR1Tails(const std::vector<std::string>& seqs,
                                                    const std::vector<std::string>& quals,
                                                    int tail_length);
    std::vector<R2HeadAnalysis> batchAnalyzeR2Heads(const std::vector<std::string>& seqs,
                                                    const std::vector<std::string>& quals,
                                                    int head_length);
    
    // 对比分析
    struct TrimComparisonResult {
        int trim_length;
        R1TailAnalysis r1_analysis;
        R2HeadAnalysis r2_analysis;
        double r1_r2_similarity_score = 0.0;  // R1尾部和R2头部的相似性评分
        double overall_quality_score = 0.0;   // 综合质量评分
    };
    
    std::vector<TrimComparisonResult> compareR1R2TrimEffects(
        const std::vector<std::string>& r1_seqs,
        const std::vector<std::string>& r2_seqs,
        const std::vector<std::string>& r1_quals,
        const std::vector<std::string>& r2_quals,
        const std::vector<int>& trim_lengths);
}

#endif