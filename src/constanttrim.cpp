#include "constanttrim.h"
#include <algorithm>
#include <iostream>
#include <cmath>

// R1尾部固定切除
int ConstantTrimmer::analyzeR1ConstantTrim(const std::string& r1_seq) {
    if (r1_seq.empty()) {
        return 0;
    }
    
    // 简单的固定切除，但需要检查序列是否足够长
    int max_possible = std::min(params.r1_trim_length, (int)r1_seq.length() - params.min_remaining_length);
    return std::max(0, max_possible);
}

// R2头部固定切除
int ConstantTrimmer::analyzeR2ConstantTrim(const std::string& r2_seq) {
    if (r2_seq.empty()) {
        return 0;
    }
    
    // 简单的固定切除，但需要检查序列是否足够长
    int max_possible = std::min(params.r2_trim_length, (int)r2_seq.length() - params.min_remaining_length);
    return std::max(0, max_possible);
}

// 基于质量的固定切除 - R1
int ConstantTrimmer::analyzeR1ConstantTrimWithQuality(const std::string& r1_seq, const std::string& r1_qual) {
    if (r1_seq.empty() || r1_qual.empty() || r1_seq.length() != r1_qual.length()) {
        return analyzeR1ConstantTrim(r1_seq);
    }
    
    int trim_length = analyzeR1ConstantTrim(r1_seq);
    if (trim_length <= 0) {
        return 0;
    }
    
    // 检查trim区域的平均质量
    if (params.use_quality_check) {
        int start_pos = r1_seq.length() - trim_length;
        double quality_sum = 0.0;
        
        for (int i = start_pos; i < (int)r1_seq.length(); i++) {
            quality_sum += (r1_qual[i] - 33); // Phred+33
        }
        
        double avg_quality = quality_sum / trim_length;
        
        // 如果trim区域质量过高，可能减少trim长度
        if (avg_quality > params.high_quality_threshold) {
            trim_length = std::max(0, trim_length - params.quality_adjustment);
        }
        // 如果trim区域质量很低，可能增加trim长度
        else if (avg_quality < params.low_quality_threshold) {
            int extra_trim = std::min(params.quality_adjustment, 
                                    (int)r1_seq.length() - start_pos - params.min_remaining_length);
            trim_length = std::min(trim_length + extra_trim, (int)r1_seq.length() - params.min_remaining_length);
        }
    }
    
    return trim_length;
}

// 基于质量的固定切除 - R2
int ConstantTrimmer::analyzeR2ConstantTrimWithQuality(const std::string& r2_seq, const std::string& r2_qual) {
    if (r2_seq.empty() || r2_qual.empty() || r2_seq.length() != r2_qual.length()) {
        return analyzeR2ConstantTrim(r2_seq);
    }
    
    int trim_length = analyzeR2ConstantTrim(r2_seq);
    if (trim_length <= 0) {
        return 0;
    }
    
    // 检查trim区域的平均质量
    if (params.use_quality_check) {
        double quality_sum = 0.0;
        
        for (int i = 0; i < trim_length; i++) {
            quality_sum += (r2_qual[i] - 33); // Phred+33
        }
        
        double avg_quality = quality_sum / trim_length;
        
        // 如果trim区域质量过高，可能减少trim长度
        if (avg_quality > params.high_quality_threshold) {
            trim_length = std::max(0, trim_length - params.quality_adjustment);
        }
        // 如果trim区域质量很低，可能增加trim长度
        else if (avg_quality < params.low_quality_threshold) {
            int extra_trim = std::min(params.quality_adjustment, 
                                    (int)r2_seq.length() - trim_length - params.min_remaining_length);
            trim_length = std::min(trim_length + extra_trim, (int)r2_seq.length() - params.min_remaining_length);
        }
    }
    
    return trim_length;
}

// 主要处理函数 - 仅序列
ConstantTrimResult ConstantTrimmer::trimReadPair(const std::string& r1_seq, const std::string& r2_seq) {
    return trimReadPair(r1_seq, r2_seq, "", "");
}

// 主要处理函数 - 序列和质量
ConstantTrimResult ConstantTrimmer::trimReadPair(const std::string& r1_seq, const std::string& r2_seq,
                                               const std::string& r1_qual, const std::string& r2_qual) {
    ConstantTrimResult result;
    
    // 分别分析R1和R2
    if (!r1_qual.empty() && r1_qual.length() == r1_seq.length()) {
        result.r1_trim_pos = analyzeR1ConstantTrimWithQuality(r1_seq, r1_qual);
    } else {
        result.r1_trim_pos = analyzeR1ConstantTrim(r1_seq);
    }
    
    if (!r2_qual.empty() && r2_qual.length() == r2_seq.length()) {
        result.r2_trim_pos = analyzeR2ConstantTrimWithQuality(r2_seq, r2_qual);
    } else {
        result.r2_trim_pos = analyzeR2ConstantTrim(r2_seq);
    }
    
    // 设置有效性
    result.r1_valid = (result.r1_trim_pos > 0);
    result.r2_valid = (result.r2_trim_pos > 0);
    
    // 计算剩余长度
    result.r1_remaining_length = r1_seq.length() - result.r1_trim_pos;
    result.r2_remaining_length = r2_seq.length() - result.r2_trim_pos;
    
    // 计算trim效率和统计
    result.r1_trim_efficiency = (double)result.r1_trim_pos / r1_seq.length() * 100.0;
    result.r2_trim_efficiency = (double)result.r2_trim_pos / r2_seq.length() * 100.0;
    
    // 分析trim区域的碱基组成
    if (result.r1_valid) {
        result.r1_trim_composition = analyzeTrimRegionComposition(r1_seq, r1_seq.length() - result.r1_trim_pos, r1_seq.length());
    }
    
    if (result.r2_valid) {
        result.r2_trim_composition = analyzeTrimRegionComposition(r2_seq, 0, result.r2_trim_pos);
    }
    
    return result;
}

// 分析trim区域碱基组成
BaseComposition ConstantTrimmer::analyzeTrimRegionComposition(const std::string& seq, int start, int end) {
    BaseComposition comp;
    
    if (start >= end || start < 0 || end > (int)seq.length()) {
        return comp;
    }
    
    int length = end - start;
    comp.total_bases = length;
    
    for (int i = start; i < end; i++) {
        char base = std::toupper(seq[i]);
        switch (base) {
            case 'A': comp.a_count++; break;
            case 'T': comp.t_count++; break;
            case 'C': comp.c_count++; break;
            case 'G': comp.g_count++; break;
            case 'N': comp.n_count++; break;
            default: comp.other_count++; break;
        }
    }
    
    // 计算百分比
    if (length > 0) {
        comp.a_percent = (double)comp.a_count / length * 100.0;
        comp.t_percent = (double)comp.t_count / length * 100.0;
        comp.c_percent = (double)comp.c_count / length * 100.0;
        comp.g_percent = (double)comp.g_count / length * 100.0;
        comp.n_percent = (double)comp.n_count / length * 100.0;
        comp.gc_percent = (double)(comp.c_count + comp.g_count) / length * 100.0;
    }
    
    return comp;
}

// 应用trim结果
std::pair<std::string, std::string> ConstantTrimmer::applyTrim(const std::string& r1_seq, 
                                                              const std::string& r2_seq, 
                                                              const ConstantTrimResult& result) {
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

// 应用trim结果到质量字符串
std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string>> 
ConstantTrimmer::applyTrimWithQuality(const std::string& r1_seq, const std::string& r2_seq,
                                     const std::string& r1_qual, const std::string& r2_qual,
                                     const ConstantTrimResult& result) {
    auto trimmed_seqs = applyTrim(r1_seq, r2_seq, result);
    
    std::string trimmed_r1_qual = r1_qual;
    std::string trimmed_r2_qual = r2_qual;
    
    // Trim质量字符串
    if (result.r1_valid && result.r1_trim_pos > 0 && result.r1_trim_pos < (int)r1_qual.length()) {
        trimmed_r1_qual = r1_qual.substr(0, r1_qual.length() - result.r1_trim_pos);
    }
    
    if (result.r2_valid && result.r2_trim_pos > 0 && result.r2_trim_pos < (int)r2_qual.length()) {
        trimmed_r2_qual = r2_qual.substr(result.r2_trim_pos);
    }
    
    return std::make_pair(trimmed_seqs, std::make_pair(trimmed_r1_qual, trimmed_r2_qual));
}

// 批量测试不同的trim长度
std::vector<ConstantTrimResult> ConstantTrimmer::testTrimLengthRange(const std::string& r1_seq, 
                                                                   const std::string& r2_seq,
                                                                   const std::string& r1_qual,
                                                                   const std::string& r2_qual,
                                                                   int min_trim, int max_trim) {
    std::vector<ConstantTrimResult> results;
    
    // 保存原始参数
    ConstantTrimParams original_params = params;
    
    for (int trim_len = min_trim; trim_len <= max_trim; trim_len++) {
        // 设置当前trim长度
        params.r1_trim_length = trim_len;
        params.r2_trim_length = trim_len;
        
        // 执行trim分析
        ConstantTrimResult result = trimReadPair(r1_seq, r2_seq, r1_qual, r2_qual);
        result.tested_trim_length = trim_len;
        
        results.push_back(result);
    }
    
    // 恢复原始参数
    params = original_params;
    
    return results;
}

// 分析不同trim长度对序列质量的影响
TrimLengthAnalysis ConstantTrimmer::analyzeTrimLengthImpact(const std::vector<std::string>& r1_seqs,
                                                           const std::vector<std::string>& r2_seqs,
                                                           const std::vector<std::string>& r1_quals,
                                                           const std::vector<std::string>& r2_quals,
                                                           int min_trim, int max_trim) {
    TrimLengthAnalysis analysis;
    
    if (r1_seqs.size() != r2_seqs.size()) {
        std::cerr << "Error: R1 and R2 sequence counts don't match!" << std::endl;
        return analysis;
    }
    
    int num_reads = r1_seqs.size();
    bool has_quality = (!r1_quals.empty() && !r2_quals.empty() && 
                       r1_quals.size() == num_reads && r2_quals.size() == num_reads);
    
    // 为每个trim长度初始化统计
    for (int trim_len = min_trim; trim_len <= max_trim; trim_len++) {
        TrimLengthStats stats;
        stats.trim_length = trim_len;
        analysis.trim_stats[trim_len] = stats;
    }
    
    // 分析每个序列对
    for (int i = 0; i < num_reads; i++) {
        std::string r1_qual = has_quality ? r1_quals[i] : "";
        std::string r2_qual = has_quality ? r2_quals[i] : "";
        
        // 测试所有trim长度
        auto results = testTrimLengthRange(r1_seqs[i], r2_seqs[i], r1_qual, r2_qual, min_trim, max_trim);
        
        for (const auto& result : results) {
            int trim_len = result.tested_trim_length;
            TrimLengthStats& stats = analysis.trim_stats[trim_len];
            
            // 累计统计
            stats.total_reads++;
            
            if (result.r1_valid) {
                stats.r1_trimmed_count++;
                stats.total_r1_bases_trimmed += result.r1_trim_pos;
                stats.total_r1_gc_in_trimmed += result.r1_trim_composition.gc_percent * result.r1_trim_pos / 100.0;
            }
            
            if (result.r2_valid) {
                stats.r2_trimmed_count++;
                stats.total_r2_bases_trimmed += result.r2_trim_pos;
                stats.total_r2_gc_in_trimmed += result.r2_trim_composition.gc_percent * result.r2_trim_pos / 100.0;
            }
            
            stats.total_remaining_r1_bases += result.r1_remaining_length;
            stats.total_remaining_r2_bases += result.r2_remaining_length;
        }
    }
    
    // 计算最终统计
    for (auto& pair : analysis.trim_stats) {
        TrimLengthStats& stats = pair.second;
        
        if (stats.total_reads > 0) {
            stats.r1_trim_rate = (double)stats.r1_trimmed_count / stats.total_reads * 100.0;
            stats.r2_trim_rate = (double)stats.r2_trimmed_count / stats.total_reads * 100.0;
            stats.avg_remaining_r1_length = (double)stats.total_remaining_r1_bases / stats.total_reads;
            stats.avg_remaining_r2_length = (double)stats.total_remaining_r2_bases / stats.total_reads;
        }
        
        if (stats.r1_trimmed_count > 0) {
            stats.avg_r1_trim_length = (double)stats.total_r1_bases_trimmed / stats.r1_trimmed_count;
            stats.avg_r1_gc_in_trimmed = stats.total_r1_gc_in_trimmed / stats.r1_trimmed_count;
        }
        
        if (stats.r2_trimmed_count > 0) {
            stats.avg_r2_trim_length = (double)stats.total_r2_bases_trimmed / stats.r2_trimmed_count;
            stats.avg_r2_gc_in_trimmed = stats.total_r2_gc_in_trimmed / stats.r2_trimmed_count;
        }
    }
    
    return analysis;
}

// 预设配置
ConstantTrimParams ConstantTrimmer::createDefaultConfig() {
    ConstantTrimParams params;
    
    params.r1_trim_length = 15;
    params.r2_trim_length = 15;
    params.min_remaining_length = 50;
    params.use_quality_check = false;
    params.high_quality_threshold = 30.0;
    params.low_quality_threshold = 10.0;
    params.quality_adjustment = 3;
    
    return params;
}

ConstantTrimParams ConstantTrimmer::createR1TailConfig(int trim_length) {
    ConstantTrimParams params = createDefaultConfig();
    
    params.r1_trim_length = trim_length;
    params.r2_trim_length = 0;  // 只trim R1
    
    return params;
}

ConstantTrimParams ConstantTrimmer::createR2HeadConfig(int trim_length) {
    ConstantTrimParams params = createDefaultConfig();
    
    params.r1_trim_length = 0;  // 只trim R2
    params.r2_trim_length = trim_length;
    
    return params;
}

ConstantTrimParams ConstantTrimmer::createQualityAwareConfig(int trim_length) {
    ConstantTrimParams params = createDefaultConfig();
    
    params.r1_trim_length = trim_length;
    params.r2_trim_length = trim_length;
    params.use_quality_check = true;
    params.high_quality_threshold = 25.0;
    params.low_quality_threshold = 15.0;
    params.quality_adjustment = 2;
    
    return params;
}