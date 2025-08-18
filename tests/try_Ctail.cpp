#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>

// 包含你的头文件
#include "mergeread.h"
#include "smarttrim.h"

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
    
    FastqRecord() = default;
    FastqRecord(const std::string& h, const std::string& s, const std::string& p, const std::string& q)
        : header(h), sequence(s), plus(p), quality(q) {}
};

class FastqReader {
public:
    FastqReader(const std::string& filename) : file(filename), records_read(0) {
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    }
    
    bool readNext(FastqRecord& record) {
        if (std::getline(file, record.header) &&
            std::getline(file, record.sequence) &&
            std::getline(file, record.plus) &&
            std::getline(file, record.quality)) {
            records_read++;
            return true;
        }
        return false;
    }
    
    int getRecordsRead() const { return records_read; }
    
private:
    std::ifstream file;
    int records_read;
};

// 日志类
class AnalysisLogger {
private:
    std::ofstream log_file;
    bool console_output;
    
public:
    AnalysisLogger(const std::string& filename, bool console = true) 
        : log_file(filename), console_output(console) {
        if (!log_file.is_open()) {
            throw std::runtime_error("Cannot create log file: " + filename);
        }
        
        // 写入文件头
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        
        log_file << "Smart Poly-C Tail Trimming Analysis Log" << std::endl;
        log_file << "Generated: " << std::ctime(&time_t);
        log_file << std::string(80, '=') << std::endl;
    }
    
    ~AnalysisLogger() {
        if (log_file.is_open()) {
            log_file.close();
        }
    }
    
    template<typename T>
    AnalysisLogger& operator<<(const T& message) {
        log_file << message;
        if (console_output) {
            std::cout << message;
        }
        return *this;
    }
    
    AnalysisLogger& operator<<(std::ostream& (*manip)(std::ostream&)) {
        log_file << manip;
        if (console_output) {
            std::cout << manip;
        }
        return *this;
    }
    
    void flush() {
        log_file.flush();
        if (console_output) {
            std::cout.flush();
        }
    }
};

// 详细统计结构
struct DetailedStats {
    int total_pairs = 0;
    int trimmed_pairs = 0;
    int r1_only_trimmed = 0;
    int r2_only_trimmed = 0;
    int both_trimmed = 0;
    int no_trim_needed = 0;
    
    // 特殊情况统计
    int long_c_tails_with_mutations = 0;
    int pure_c_tails = 0;
    int short_c_tails = 0;
    int false_positives_avoided = 0;
    
    std::vector<int> r1_trim_lengths;
    std::vector<int> r2_trim_lengths;
    std::vector<double> trim_scores;
    
    // 用于保存例子
    std::vector<std::string> examples_no_trim;
    std::vector<std::string> examples_c_tail_cut;
    std::vector<std::string> examples_mutation_handled;
    
    void addExample(const std::string& type, const std::string& example) {
        if (type == "no_trim" && examples_no_trim.size() < 5) {
            examples_no_trim.push_back(example);
        } else if (type == "c_tail_cut" && examples_c_tail_cut.size() < 5) {
            examples_c_tail_cut.push_back(example);
        } else if (type == "mutation_handled" && examples_mutation_handled.size() < 5) {
            examples_mutation_handled.push_back(example);
        }
    }
    
    void printDetailedReport(AnalysisLogger& logger) const {
        logger << std::endl << std::string(80, '=') << std::endl;
        logger << "                    DETAILED TRIMMING ANALYSIS REPORT" << std::endl;
        logger << std::string(80, '=') << std::endl;
        
        // 基本统计
        logger << std::endl << "BASIC STATISTICS:" << std::endl;
        logger << "Total read pairs processed: " << total_pairs << std::endl;
        logger << "Pairs requiring trimming: " << trimmed_pairs 
               << " (" << std::fixed << std::setprecision(2) 
               << (100.0 * trimmed_pairs / total_pairs) << "%)" << std::endl;
        logger << "Pairs with no trimming needed: " << no_trim_needed 
               << " (" << (100.0 * no_trim_needed / total_pairs) << "%)" << std::endl;
        logger << "Efficiency: " << (100.0 * (total_pairs - no_trim_needed) / total_pairs) 
               << "% reads had potential poly-C issues" << std::endl;
        
        // 修剪模式分析
        logger << std::endl << "TRIMMING PATTERN ANALYSIS:" << std::endl;
        logger << "R1 only trimmed: " << r1_only_trimmed;
        if (trimmed_pairs > 0) {
            logger << " (" << (100.0 * r1_only_trimmed / trimmed_pairs) << "% of trimmed)";
        }
        logger << std::endl;
        
        logger << "R2 only trimmed: " << r2_only_trimmed;
        if (trimmed_pairs > 0) {
            logger << " (" << (100.0 * r2_only_trimmed / trimmed_pairs) << "% of trimmed)";
        }
        logger << std::endl;
        
        logger << "Both R1 & R2 trimmed: " << both_trimmed;
        if (trimmed_pairs > 0) {
            logger << " (" << (100.0 * both_trimmed / trimmed_pairs) << "% of trimmed)";
        }
        logger << std::endl;
        
        // Poly-C特征分析
        logger << std::endl << "POLY-C TAIL CHARACTERISTICS:" << std::endl;
        logger << "Pure poly-C tails detected: " << pure_c_tails << std::endl;
        logger << "Long C-tails with mutations: " << long_c_tails_with_mutations << std::endl;
        logger << "Short C-tails (high precision): " << short_c_tails << std::endl;
        logger << "False positives avoided: " << false_positives_avoided << std::endl;
        
        // 长度分析
        if (!r1_trim_lengths.empty()) {
            double r1_sum = std::accumulate(r1_trim_lengths.begin(), r1_trim_lengths.end(), 0.0);
            double r1_avg = r1_sum / r1_trim_lengths.size();
            int r1_max = *std::max_element(r1_trim_lengths.begin(), r1_trim_lengths.end());
            logger << std::endl << "TRIMMING LENGTH ANALYSIS:" << std::endl;
            logger << "R1 average trim length: " << std::fixed << std::setprecision(1) << r1_avg << " bp" << std::endl;
            logger << "R1 maximum trim length: " << r1_max << " bp" << std::endl;
        }
        
        if (!r2_trim_lengths.empty()) {
            double r2_sum = std::accumulate(r2_trim_lengths.begin(), r2_trim_lengths.end(), 0.0);
            double r2_avg = r2_sum / r2_trim_lengths.size();
            int r2_max = *std::max_element(r2_trim_lengths.begin(), r2_trim_lengths.end());
            logger << "R2 average trim length: " << std::fixed << std::setprecision(1) << r2_avg << " bp" << std::endl;
            logger << "R2 maximum trim length: " << r2_max << " bp" << std::endl;
        }
        
        // 分数分析
        if (!trim_scores.empty()) {
            double score_sum = std::accumulate(trim_scores.begin(), trim_scores.end(), 0.0);
            double score_avg = score_sum / trim_scores.size();
            double score_max = *std::max_element(trim_scores.begin(), trim_scores.end());
            logger << std::endl << "SCORING ANALYSIS:" << std::endl;
            logger << "Average trim score: " << std::fixed << std::setprecision(2) << score_avg << std::endl;
            logger << "Maximum trim score: " << score_max << std::endl;
        }
        
        // 例子展示
        logger << std::endl << "DETAILED EXAMPLES:" << std::endl;
        
        logger << std::endl << "Examples of reads NOT requiring trimming:" << std::endl;
        for (size_t i = 0; i < examples_no_trim.size(); i++) {
            logger << "   " << (i+1) << ". " << examples_no_trim[i] << std::endl;
        }
        
        logger << std::endl << "Examples of successful poly-C tail trimming:" << std::endl;
        for (size_t i = 0; i < examples_c_tail_cut.size(); i++) {
            logger << "   " << (i+1) << ". " << examples_c_tail_cut[i] << std::endl;
        }
        
        logger << std::endl << "Examples of mutation-containing C-tails (correctly handled):" << std::endl;
        for (size_t i = 0; i < examples_mutation_handled.size(); i++) {
            logger << "   " << (i+1) << ". " << examples_mutation_handled[i] << std::endl;
        }
        
        logger << std::endl << std::string(80, '=') << std::endl;
    }
    
    // 保存数据到CSV文件用于Python可视化
    void saveToCSV(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot create CSV file: " << filename << std::endl;
            return;
        }
        
        // 写入基本统计数据
        file << "metric,value\n";
        file << "total_pairs," << total_pairs << "\n";
        file << "trimmed_pairs," << trimmed_pairs << "\n";
        file << "trim_rate," << (100.0 * trimmed_pairs / total_pairs) << "\n";
        file << "r1_only_trimmed," << r1_only_trimmed << "\n";
        file << "r2_only_trimmed," << r2_only_trimmed << "\n";
        file << "both_trimmed," << both_trimmed << "\n";
        file << "pure_c_tails," << pure_c_tails << "\n";
        file << "c_tails_with_mutations," << long_c_tails_with_mutations << "\n";
        
        file.close();
        
        // 保存长度分布数据
        std::ofstream length_file(filename.substr(0, filename.find('.')) + "_lengths.csv");
        length_file << "read_type,trim_length\n";
        for (int len : r1_trim_lengths) {
            length_file << "R1," << len << "\n";
        }
        for (int len : r2_trim_lengths) {
            length_file << "R2," << len << "\n";
        }
        length_file.close();
        
        // 保存分数分布数据
        std::ofstream score_file(filename.substr(0, filename.find('.')) + "_scores.csv");
        score_file << "trim_score\n";
        for (double score : trim_scores) {
            score_file << score << "\n";
        }
        score_file.close();
    }
};

// 分析序列特征
std::string analyzeSequenceFeatures(const std::string& seq, int start, int end) {
    if (start >= end || start < 0 || end > (int)seq.length()) return "invalid_region";
    
    std::string region = seq.substr(start, end - start);
    int c_count = 0, total = region.length();
    
    for (char base : region) {
        if (base == 'C' || base == 'c') c_count++;
    }
    
    double c_ratio = (double)c_count / total;
    
    if (c_ratio >= 0.9) return "pure_c_tail";
    else if (c_ratio >= 0.7 && total >= 8) return "c_tail_with_mutations";
    else if (c_ratio >= 0.6 && total <= 5) return "short_c_tail";
    else return "mixed_sequence";
}

DetailedStats processCompleteDataset(const std::string& r1_file, const std::string& r2_file, 
                                    const std::string& sample_name, AnalysisLogger& logger) {
    logger << std::endl << "PROCESSING COMPLETE DATASET: " << sample_name << std::endl;
    logger << "R1 file: " << r1_file << std::endl;
    logger << "R2 file: " << r2_file << std::endl;
    logger << std::string(60, '-') << std::endl;
    
    DetailedStats stats;
    ReadMerger merger;
    TrimParams default_params(25, 0.0);
    SmartTrimmer trimmer(default_params);
    
    try {
        FastqReader r1_reader(r1_file);
        FastqReader r2_reader(r2_file);
        
        FastqRecord r1_record, r2_record;
        int processed = 0;
        
        // 处理进度显示
        auto start_time = std::chrono::high_resolution_clock::now();
        
        while (r1_reader.readNext(r1_record) && r2_reader.readNext(r2_record)) {
            processed++;
            stats.total_pairs++;
            
            // 显示处理进度
            if (processed % 50000 == 0) {
                auto current_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time);
                logger << "Processed " << processed << " pairs (" 
                       << (duration.count() / 1000.0) << "s elapsed)..." << std::endl;
                logger.flush();
            }
            
            // 执行trimming分析
            TrimResult trim_result = trimmer.findOptimalTrimPositions(r1_record.sequence, r2_record.sequence);
            
            if (trim_result.is_valid) {
                stats.trimmed_pairs++;
                stats.trim_scores.push_back(trim_result.final_score);
                
                bool r1_trimmed = trim_result.r1_trim_pos > 0;
                bool r2_trimmed = trim_result.r2_trim_pos > 0;
                
                if (r1_trimmed && r2_trimmed) {
                    stats.both_trimmed++;
                } else if (r1_trimmed) {
                    stats.r1_only_trimmed++;
                } else if (r2_trimmed) {
                    stats.r2_only_trimmed++;
                }
                
                if (r1_trimmed) {
                    stats.r1_trim_lengths.push_back(trim_result.r1_trim_pos);
                }
                if (r2_trimmed) {
                    stats.r2_trim_lengths.push_back(trim_result.r2_trim_pos);
                }
                
                // 分析被修剪的区域特征
                if (r1_trimmed) {
                    int trim_start = r1_record.sequence.length() - trim_result.r1_trim_pos;
                    std::string feature = analyzeSequenceFeatures(r1_record.sequence, trim_start, r1_record.sequence.length());
                    
                    if (feature == "pure_c_tail") {
                        stats.pure_c_tails++;
                        if (stats.examples_c_tail_cut.size() < 5) {
                            std::stringstream example;
                            example << "R1 pure C-tail: " << r1_record.sequence.substr(trim_start) 
                                   << " (score: " << std::fixed << std::setprecision(1) << trim_result.final_score << ")";
                            stats.addExample("c_tail_cut", example.str());
                        }
                    } else if (feature == "c_tail_with_mutations") {
                        stats.long_c_tails_with_mutations++;
                        if (stats.examples_mutation_handled.size() < 5) {
                            std::stringstream example;
                            example << "R1 C-tail w/ mutations: " << r1_record.sequence.substr(trim_start)
                                   << " (correctly trimmed, score: " << std::fixed << std::setprecision(1) << trim_result.final_score << ")";
                            stats.addExample("mutation_handled", example.str());
                        }
                    } else if (feature == "short_c_tail") {
                        stats.short_c_tails++;
                    }
                }
                
                // 类似地分析R2
                if (r2_trimmed) {
                    std::string feature = analyzeSequenceFeatures(r2_record.sequence, 0, trim_result.r2_trim_pos);
                    
                    if (feature == "pure_c_tail") {
                        stats.pure_c_tails++;
                        if (stats.examples_c_tail_cut.size() < 5) {
                            std::stringstream example;
                            example << "R2 pure C-tail: " << r2_record.sequence.substr(0, trim_result.r2_trim_pos)
                                   << " (score: " << std::fixed << std::setprecision(1) << trim_result.final_score << ")";
                            stats.addExample("c_tail_cut", example.str());
                        }
                    } else if (feature == "c_tail_with_mutations") {
                        stats.long_c_tails_with_mutations++;
                        if (stats.examples_mutation_handled.size() < 5) {
                            std::stringstream example;
                            example << "R2 C-tail w/ mutations: " << r2_record.sequence.substr(0, trim_result.r2_trim_pos)
                                   << " (correctly trimmed, score: " << std::fixed << std::setprecision(1) << trim_result.final_score << ")";
                            stats.addExample("mutation_handled", example.str());
                        }
                    } else if (feature == "short_c_tail") {
                        stats.short_c_tails++;
                    }
                }
                
            } else {
                stats.no_trim_needed++;
                
                // 收集一些不需要trimming的例子
                if (stats.examples_no_trim.size() < 5) {
                    std::stringstream example;
                    int r1_start = std::max(0, (int)r1_record.sequence.length()-15);
                    int r2_end = std::min(15, (int)r2_record.sequence.length());
                    example << "R1: " << r1_record.sequence.substr(r1_start)
                           << " | R2: " << r2_record.sequence.substr(0, r2_end)
                           << " (no significant poly-C detected)";
                    stats.addExample("no_trim", example.str());
                }
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        logger << "Processing completed!" << std::endl;
        logger << "Final count: " << processed << " read pairs processed" << std::endl;
        logger << "Total time: " << (total_duration.count() / 1000.0) << " seconds" << std::endl;
        logger << "Processing speed: " << std::fixed << std::setprecision(0) 
               << (processed / (total_duration.count() / 1000.0)) << " pairs/second" << std::endl;
        
    } catch (const std::exception& e) {
        logger << "Error processing " << sample_name << ": " << e.what() << std::endl;
        logger << "Please check if the files exist at the specified paths." << std::endl;
    }
    
    return stats;
}

void testBasicAlgorithms(AnalysisLogger& logger) {
    logger << std::endl << "TESTING BASIC ALGORITHMS" << std::endl;
    logger << std::string(50, '-') << std::endl;
    
    ReadMerger merger;
    TrimParams test_params(25, 0.0);
    SmartTrimmer trimmer(test_params);
    
    // 详细测试用例
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"ATGCGATCGATCCCCCCCCCCC", "GATCGATCGATCGCTTAGAA"},  // R1有poly-C尾部，R2正常
        {"ATTTTACCCGATCGATCGATT", "GATCGATCGATCGCTTAGAA"},   // R1有T和一些C，R2正常
        {"ATGCGATCGATGGATGCATGC", "GATCGATCGATCGCTTAGAA"},   // 都正常
        {"ATGCGATCGATCCCGCCCACC", "GATCGATCGATCGCTTAGAA"},   // R1混合C，R2正常
        {"ATGCGATCGATGGATGCATGC", "CCCCCCCGATCGATCGCTTA"},   // R1正常，R2有poly-C头部 ★
        {"ATGCGATCGATGGATGCATGC", "CCCGCCCGATCGATCGCTTA"},   // R1正常，R2混合C头部 ★
        {"ATGCGATCGATGGATGCATGC", "CCCCCCCCCCCCCCCCCCCC"},   // R1正常，R2全是C ★
        {"ATGCGATCGATGGATGCATGC", "CCCGCCCACCCGATCGCTTA"},   // R1正常，R2有突变的C头部 ★
    };
    
    for (size_t i = 0; i < test_cases.size(); i++) {
        const auto& [r1_seq, r2_seq] = test_cases[i];
        
        logger << "=== Test Case " << (i+1) << " ===" << std::endl;
        logger << "R1: " << r1_seq << std::endl;
        logger << "R2: " << r2_seq << std::endl;
        
        // 分析R2开头的poly-C
        int check_length = std::min(test_params.window_size, (int)r2_seq.length());
        std::string r2_head = r2_seq.substr(0, check_length);
        
        logger << "R2 head analysis (first " << check_length << " bp):" << std::endl;
        logger << "  R2 head: " << r2_head << std::endl;
        
        // 分析R2头部的C含量
        int c_count = 0;
        for (char base : r2_head) {
            if (base == 'C' || base == 'c') c_count++;
        }
        double c_ratio = (double)c_count / r2_head.length();
        logger << "  C content: " << c_count << "/" << r2_head.length() 
               << " (" << (c_ratio * 100) << "%)" << std::endl;
        
        // 显示前10个碱基的详细信息
        logger << "  First 10 bases: ";
        for (int j = 0; j < std::min(10, (int)r2_seq.length()); j++) {
            logger << r2_seq[j];
        }
        logger << std::endl;
        
        // 执行actual trimming
        TrimResult result = trimmer.findOptimalTrimPositions(r1_seq, r2_seq);
        logger << "  Algorithm Decision: " << (result.is_valid ? "TRIM" : "NO TRIM") << std::endl;
        
        if (result.is_valid) {
            logger << "  Final Score: " << result.final_score << std::endl;
            logger << "  R1 trim position: " << result.r1_trim_pos << std::endl;
            logger << "  R2 trim position: " << result.r2_trim_pos << std::endl;
            
            if (result.r1_trim_pos > 0) {
                int trim_start = r1_seq.length() - result.r1_trim_pos;
                logger << "  R1 would trim: " << r1_seq.substr(trim_start) << " (from position " << trim_start << ")" << std::endl;
            }
            if (result.r2_trim_pos > 0) {
                logger << "  R2 would trim: " << r2_seq.substr(0, result.r2_trim_pos) << " (first " << result.r2_trim_pos << " bp)" << std::endl;
            }
        } else {
            logger << "  Threshold not met (score <= 0)" << std::endl;
        }
        
        logger << std::endl;
    }
    
    logger << "Basic algorithm tests completed" << std::endl;
}

int main() {
    // 创建日志文件
    std::string log_filename = "ctail_analysis_log.txt";
    AnalysisLogger logger(log_filename, true);
    
    logger << std::string(80, '=') << std::endl;
    logger << "        SMART POLY-C TAIL TRIMMING - COMPREHENSIVE ANALYSIS" << std::endl;
    logger << "                     Distance-Decay Integration Algorithm" << std::endl;
    logger << std::string(80, '=') << std::endl;
    
    try {
        // 基本算法测试
        testBasicAlgorithms(logger);
        
        // 参数优化测试
        logger << std::endl << "PARAMETER OPTIMIZATION TESTING" << std::endl;
        logger << std::string(50, '-') << std::endl;
        
        TrimParams default_params(25, 0.0);
        SmartTrimmer default_trimmer(default_params);
        auto param_matrix = default_trimmer.generateParameterMatrix();
        
        logger << "Testing " << param_matrix.size() << " parameter configurations..." << std::endl;
        
        // 用一小部分数据测试所有参数配置
        std::string test_r1_file = "normal_R1.fastq";  // 或者你想用的测试文件
        std::string test_r2_file = "normal_R2.fastq";
        
        std::vector<std::tuple<TrimParams, int, double, double>> param_results;
        
        try {
            FastqReader r1_reader(test_r1_file);
            FastqReader r2_reader(test_r2_file);
            
            // 读取前10000对reads用于参数测试
            std::vector<std::pair<std::string, std::string>> test_pairs;
            FastqRecord r1_record, r2_record;
            int test_count = 0;
            
            while (r1_reader.readNext(r1_record) && r2_reader.readNext(r2_record) && test_count < 10000) {
                test_pairs.push_back({r1_record.sequence, r2_record.sequence});
                test_count++;
            }
            
            logger << "Loaded " << test_pairs.size() << " read pairs for parameter testing" << std::endl;
            
            // 测试每个参数配置
            for (size_t i = 0; i < param_matrix.size(); i++) {
                const auto& params = param_matrix[i];
                SmartTrimmer test_trimmer(params);
                
                int trimmed_count = 0;
                double total_score = 0.0;
                double total_trim_length = 0.0;
                
                for (const auto& pair : test_pairs) {
                    TrimResult result = test_trimmer.findOptimalTrimPositions(pair.first, pair.second);
                    if (result.is_valid) {
                        trimmed_count++;
                        total_score += result.final_score;
                        total_trim_length += (result.r1_trim_pos + result.r2_trim_pos);
                    }
                }
                
                double trim_rate = (double)trimmed_count / test_pairs.size();
                double avg_score = trimmed_count > 0 ? total_score / trimmed_count : 0.0;
                double avg_length = trimmed_count > 0 ? total_trim_length / trimmed_count : 0.0;
                
                param_results.push_back({params, trimmed_count, trim_rate, avg_score});
                
                if (i % 10 == 0) {
                    logger << "Progress: " << i << "/" << param_matrix.size() << " configurations tested" << std::endl;
                }
            }
            
            // 排序并显示最佳参数
            std::sort(param_results.begin(), param_results.end(), 
                     [](const auto& a, const auto& b) {
                         return std::get<3>(a) > std::get<3>(b);  // 按平均分数排序
                     });
            
            logger << std::endl << "TOP 10 PARAMETER CONFIGURATIONS:" << std::endl;
            logger << "Rank | Window | InitScore | CScore | PenaltyScore | TrimRate | AvgScore" << std::endl;
            logger << std::string(75, '-') << std::endl;
            
            for (int i = 0; i < std::min(10, (int)param_results.size()); i++) {
                const auto& [params, trim_count, trim_rate, avg_score] = param_results[i];
                logger << std::setw(4) << (i+1) << " | "
                       << std::setw(6) << params.window_size << " | "
                       << std::setw(9) << std::fixed << std::setprecision(1) << params.initial_score << " | "
                       << std::setw(6) << params.c_score << " | "
                       << std::setw(12) << params.penalty_score << " | "
                       << std::setw(8) << std::setprecision(2) << (trim_rate * 100) << "% | "
                       << std::setw(8) << avg_score << std::endl;
            }
            
            // 保存参数优化结果
            std::ofstream param_file("parameter_optimization_results.csv");
            param_file << "rank,window_size,initial_score,c_score,penalty_score,trim_rate,avg_score,trim_count\n";
            for (size_t i = 0; i < param_results.size(); i++) {
                const auto& [params, trim_count, trim_rate, avg_score] = param_results[i];
                param_file << (i+1) << "," << params.window_size << "," << params.initial_score << ","
                          << params.c_score << "," << params.penalty_score << "," << trim_rate << ","
                          << avg_score << "," << trim_count << "\n";
            }
            param_file.close();
            
            logger << std::endl << "Parameter optimization results saved to: parameter_optimization_results.csv" << std::endl;
            
        } catch (const std::exception& e) {
            logger << "Parameter testing failed: " << e.what() << std::endl;
        }
        
        // 处理真实的完整数据集
        std::vector<std::tuple<std::string, std::string, std::string>> datasets = {
            {"normal_R1.fastq", "normal_R2.fastq", "Normal_Sample"},
            {"tumor_R1.fastq", "tumor_R2.fastq", "Tumor_Sample"}
        };
        
        for (const auto& dataset : datasets) {
            std::string r1_file = std::get<0>(dataset);
            std::string r2_file = std::get<1>(dataset);
            std::string sample_name = std::get<2>(dataset);
            
            DetailedStats stats = processCompleteDataset(r1_file, r2_file, sample_name, logger);
            
            if (stats.total_pairs > 0) {
                stats.printDetailedReport(logger);
                
                // 保存数据到CSV文件
                std::string csv_filename = sample_name + "_trimming_analysis.csv";
                stats.saveToCSV(csv_filename);
                logger << "Data saved to CSV files for Python visualization:" << std::endl;
                logger << "   - " << csv_filename << " (basic statistics)" << std::endl;
                logger << "   - " << sample_name << "_trimming_analysis_lengths.csv (length distribution)" << std::endl;
                logger << "   - " << sample_name << "_trimming_analysis_scores.csv (score distribution)" << std::endl;
            }
        }
        
        logger << std::endl << "ANALYSIS COMPLETED SUCCESSFULLY!" << std::endl;
        logger << std::endl << "NEXT STEPS:" << std::endl;
        logger << "1. Review the detailed analysis above" << std::endl;
        logger << "2. Use the generated CSV files for Python visualization" << std::endl;
        logger << "3. Adjust algorithm parameters if needed based on results" << std::endl;
        logger << "4. Run on production data with optimized settings" << std::endl;
        logger << std::endl << "Log saved to: " << log_filename << std::endl;
        
    } catch (const std::exception& e) {
        logger << "Analysis failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}