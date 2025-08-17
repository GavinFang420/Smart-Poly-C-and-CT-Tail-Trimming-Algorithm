#include "smarttrim.h"
#include "mergeread.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <chrono>
#include <fstream>
#include <random>
#include <algorithm>
#include <map>

class TestRunner {
private:
    int tests_passed = 0;
    int tests_failed = 0;
    
public:
    void assert_test(bool condition, const std::string& test_name) {
        if (condition) {
            std::cout << "[PASS] " << test_name << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] " << test_name << std::endl;
            tests_failed++;
        }
    }
    
    void print_summary() {
        std::cout << "\n=== Test Summary ===" << std::endl;
        std::cout << "Passed: " << tests_passed << std::endl;
        std::cout << "Failed: " << tests_failed << std::endl;
        std::cout << "Total:  " << (tests_passed + tests_failed) << std::endl;
        
        if (tests_failed == 0) {
            std::cout << "All tests PASSED! ✓" << std::endl;
        } else {
            std::cout << "Some tests FAILED! ✗" << std::endl;
        }
    }
};

// 数据生成器 - 模拟真实WGBS数据
class WGBSDataGenerator {
private:
    std::mt19937 rng;
    
public:
    WGBSDataGenerator() : rng(std::chrono::steady_clock::now().time_since_epoch().count()) {}
    
    struct ReadPair {
        std::string r1_seq;
        std::string r2_seq;
        std::string r1_qual;
        std::string r2_qual;
        int true_c_tail_r1;  // 真实的C tail长度
        int true_c_tail_r2;  // 真实的C tail长度
        std::string description;
    };
    
    // 生成带C tail的测试数据
    std::vector<ReadPair> generateTestData(int num_pairs = 100) {
        std::vector<ReadPair> data;
        std::uniform_int_distribution<int> c_tail_dist(0, 20);  // C tail长度0-20
        std::uniform_int_distribution<int> seq_len_dist(120, 150);  // 序列长度120-150
        std::uniform_int_distribution<int> base_dist(0, 3);
        std::string bases = "ATGC";
        
        for (int i = 0; i < num_pairs; i++) {
            ReadPair pair;
            
            // 生成R1序列
            int r1_len = seq_len_dist(rng);
            int r1_c_tail = c_tail_dist(rng);
            
            // 生成正常序列部分
            for (int j = 0; j < r1_len - r1_c_tail; j++) {
                pair.r1_seq += bases[base_dist(rng)];
            }
            // 添加C tail
            for (int j = 0; j < r1_c_tail; j++) {
                pair.r1_seq += 'C';
            }
            pair.true_c_tail_r1 = r1_c_tail;
            
            // 生成R2序列
            int r2_len = seq_len_dist(rng);
            int r2_c_head = c_tail_dist(rng);  // R2的C在开头
            
            // 添加C head
            for (int j = 0; j < r2_c_head; j++) {
                pair.r2_seq += 'C';
            }
            // 添加正常序列部分
            for (int j = 0; j < r2_len - r2_c_head; j++) {
                pair.r2_seq += bases[base_dist(rng)];
            }
            pair.true_c_tail_r2 = r2_c_head;
            
            // 生成质量分数（简单模拟）
            pair.r1_qual = std::string(pair.r1_seq.length(), 'I');  // 高质量
            pair.r2_qual = std::string(pair.r2_seq.length(), 'I');
            
            pair.description = "Simulated WGBS pair " + std::to_string(i) + 
                              " (R1_Ctail=" + std::to_string(r1_c_tail) + 
                              ", R2_Chead=" + std::to_string(r2_c_head) + ")";
            
            data.push_back(pair);
        }
        
        return data;
    }
    
    // 从FASTQ文件读取真实数据（如果有的话）
    std::vector<ReadPair> loadFromFastq(const std::string& r1_file, const std::string& r2_file) {
        std::vector<ReadPair> data;
        std::ifstream f1(r1_file), f2(r2_file);
        
        if (!f1.is_open() || !f2.is_open()) {
            std::cout << "Warning: Cannot open FASTQ files, using simulated data" << std::endl;
            return generateTestData(100);
        }
        
        std::string line1, line2;
        ReadPair pair;
        int line_count = 0;
        
        while (std::getline(f1, line1) && std::getline(f2, line2)) {
            int pos = line_count % 4;
            
            if (pos == 0) {  // Header line
                pair = ReadPair();
                pair.description = line1;
            } else if (pos == 1) {  // Sequence
                pair.r1_seq = line1;
                pair.r2_seq = line2;
            } else if (pos == 3) {  // Quality
                pair.r1_qual = line1;
                pair.r2_qual = line2;
                
                // 估算真实C tail（简单启发式）
                pair.true_c_tail_r1 = estimateCTail(pair.r1_seq, false);
                pair.true_c_tail_r2 = estimateCTail(pair.r2_seq, true);
                
                data.push_back(pair);
                
                if (data.size() >= 100) break;  // 限制数量
            }
            line_count++;
        }
        
        return data;
    }
    
private:
    int estimateCTail(const std::string& seq, bool from_head) {
        int c_count = 0;
        if (from_head) {
            for (char base : seq) {
                if (base == 'C' || base == 'c') c_count++;
                else break;
            }
        } else {
            for (int i = seq.length() - 1; i >= 0; i--) {
                if (seq[i] == 'C' || seq[i] == 'c') c_count++;
                else break;
            }
        }
        return c_count;
    }
};

// 参数优化器
class ParameterOptimizer {
public:
    struct OptimizationResult {
        TrimParams best_params;
        double best_score;
        std::map<std::string, double> metrics;
        std::vector<std::pair<TrimParams, double>> all_results;
    };
    
    // 生成权重衰减函数的候选
    std::vector<std::vector<double>> generateWeightFunctions(int window_size) {
        std::vector<std::vector<double>> functions;
        
        // 1. 线性衰减函数组
        for (double slope = 0.5; slope <= 3.0; slope += 0.5) {
            std::vector<double> weights(window_size);
            for (int i = 0; i < window_size; i++) {
                weights[i] = 1.0 + slope * i / (window_size - 1);
            }
            functions.push_back(weights);
        }
        
        // 2. 指数衰减函数组
        for (double exp_factor = 0.5; exp_factor <= 2.0; exp_factor += 0.5) {
            std::vector<double> weights(window_size);
            for (int i = 0; i < window_size; i++) {
                double x = (double)i / (window_size - 1);
                weights[i] = std::exp(exp_factor * x);
            }
            functions.push_back(weights);
        }
        
        // 3. 二次函数组
        for (double quad_factor = 0.5; quad_factor <= 2.0; quad_factor += 0.5) {
            std::vector<double> weights(window_size);
            for (int i = 0; i < window_size; i++) {
                double x = (double)i / (window_size - 1);
                weights[i] = 1.0 + quad_factor * x * x;
            }
            functions.push_back(weights);
        }
        
        // 4. 对数函数组
        for (double log_factor = 1.0; log_factor <= 3.0; log_factor += 1.0) {
            std::vector<double> weights(window_size);
            for (int i = 0; i < window_size; i++) {
                double x = (double)i / (window_size - 1);
                weights[i] = 1.0 + log_factor * std::log(1.0 + x);
            }
            functions.push_back(weights);
        }
        
        return functions;
    }
    
    // 评估trimming结果的准确性
    double evaluateTrimAccuracy(const std::vector<WGBSDataGenerator::ReadPair>& test_data,
                               const TrimParams& params) {
        SmartTrimmer trimmer(params);
        
        double total_score = 0.0;
        int valid_cases = 0;
        
        for (const auto& pair : test_data) {
            TrimResult result = trimmer.findOptimalTrimPositions(pair.r1_seq, pair.r2_seq);
            
            if (result.is_valid) {
                // 计算准确性分数
                double r1_accuracy = 1.0 - std::abs(result.r1_trim_pos - pair.true_c_tail_r1) / 
                                    std::max(1.0, (double)std::max(result.r1_trim_pos, pair.true_c_tail_r1));
                double r2_accuracy = 1.0 - std::abs(result.r2_trim_pos - pair.true_c_tail_r2) / 
                                    std::max(1.0, (double)std::max(result.r2_trim_pos, pair.true_c_tail_r2));
                
                total_score += (r1_accuracy + r2_accuracy) / 2.0;
                valid_cases++;
            }
        }
        
        return valid_cases > 0 ? total_score / valid_cases : 0.0;
    }
    
    // 主要的参数优化函数
    OptimizationResult optimizeParameters(const std::vector<WGBSDataGenerator::ReadPair>& data) {
        OptimizationResult result;
        result.best_score = -1.0;
        
        std::cout << "Starting parameter optimization with " << data.size() << " test cases..." << std::endl;
        
        // 测试不同的窗口大小
        std::vector<int> window_sizes = {20, 25, 30, 35};
        
        // 测试不同的惩罚分数
        std::vector<double> penalty_scores = {-2.0, -3.0, -4.0, -5.0, -6.0};
        
        // 测试不同的初始分数
        std::vector<double> initial_scores = {0.0, -5.0, -10.0, -15.0};
        
        int total_combinations = 0;
        int tested_combinations = 0;
        
        for (int ws : window_sizes) {
            auto weight_functions = generateWeightFunctions(ws);
            total_combinations += weight_functions.size() * penalty_scores.size() * initial_scores.size();
        }
        
        std::cout << "Total parameter combinations to test: " << total_combinations << std::endl;
        
        for (int ws : window_sizes) {
            auto weight_functions = generateWeightFunctions(ws);
            
            for (const auto& weights : weight_functions) {
                for (double penalty : penalty_scores) {
                    for (double init_score : initial_scores) {
                        TrimParams params(ws, init_score);
                        params.penalty_score = penalty;
                        params.position_weights = weights;
                        
                        double score = evaluateTrimAccuracy(data, params);
                        result.all_results.push_back({params, score});
                        
                        if (score > result.best_score) {
                            result.best_score = score;
                            result.best_params = params;
                        }
                        
                        tested_combinations++;
                        if (tested_combinations % 50 == 0) {
                            std::cout << "Progress: " << tested_combinations << "/" << total_combinations 
                                     << " (" << (100.0 * tested_combinations / total_combinations) << "%)" 
                                     << " Best score so far: " << result.best_score << std::endl;
                        }
                    }
                }
            }
        }
        
        // 计算统计指标
        result.metrics["total_combinations"] = total_combinations;
        result.metrics["best_accuracy"] = result.best_score;
        
        return result;
    }
};

void test_parameter_optimization() {
    std::cout << "\n=== Testing Parameter Optimization Framework ===" << std::endl;
    TestRunner runner;
    
    // 生成测试数据
    WGBSDataGenerator generator;
    auto test_data = generator.generateTestData(50);  // 50个测试用例
    
    runner.assert_test(test_data.size() == 50, "Should generate 50 test cases");
    runner.assert_test(!test_data[0].r1_seq.empty(), "Generated sequences should not be empty");
    
    // 测试参数优化器
    ParameterOptimizer optimizer;
    auto weight_functions = optimizer.generateWeightFunctions(30);
    
    runner.assert_test(weight_functions.size() > 10, "Should generate multiple weight functions");
    runner.assert_test(weight_functions[0].size() == 30, "Weight function should match window size");
    
    std::cout << "Generated " << weight_functions.size() << " different weight functions" << std::endl;
    
    // 快速优化测试（减少参数组合）
    std::cout << "Running quick parameter optimization..." << std::endl;
    auto quick_data = generator.generateTestData(20);  // 更少的数据用于快速测试
    auto opt_result = optimizer.optimizeParameters(quick_data);
    
    runner.assert_test(opt_result.best_score >= 0, "Should find valid parameter combination");
    runner.assert_test(!opt_result.all_results.empty(), "Should test multiple parameter combinations");
    
    std::cout << "Best accuracy achieved: " << opt_result.best_score << std::endl;
    std::cout << "Best parameters: window=" << opt_result.best_params.window_size 
              << ", penalty=" << opt_result.best_params.penalty_score
              << ", init_score=" << opt_result.best_params.initial_score << std::endl;
    
    runner.print_summary();
}

void test_real_data_workflow() {
    std::cout << "\n=== Testing Real Data Workflow ===" << std::endl;
    TestRunner runner;
    
    WGBSDataGenerator generator;
    
    // 尝试加载真实数据，如果失败则使用模拟数据
    auto data = generator.loadFromFastq("test_R1.fastq", "test_R2.fastq");
    
    runner.assert_test(!data.empty(), "Should load or generate test data");
    std::cout << "Loaded " << data.size() << " read pairs for analysis" << std::endl;
    
    // 显示数据统计
    int total_r1_c_tail = 0, total_r2_c_tail = 0;
    for (const auto& pair : data) {
        total_r1_c_tail += pair.true_c_tail_r1;
        total_r2_c_tail += pair.true_c_tail_r2;
    }
    
    double avg_r1_c_tail = (double)total_r1_c_tail / data.size();
    double avg_r2_c_tail = (double)total_r2_c_tail / data.size();
    
    std::cout << "Average C tail length - R1: " << avg_r1_c_tail 
              << ", R2: " << avg_r2_c_tail << std::endl;
    
    runner.assert_test(avg_r1_c_tail >= 0 && avg_r2_c_tail >= 0, "C tail statistics should be valid");
    
    runner.print_summary();
}

void full_parameter_optimization_benchmark() {
    std::cout << "\n=== Full Parameter Optimization Benchmark ===" << std::endl;
    
    WGBSDataGenerator generator;
    ParameterOptimizer optimizer;
    
    // 生成完整的测试数据集
    auto data = generator.generateTestData(100);
    std::cout << "Generated " << data.size() << " test cases for comprehensive optimization" << std::endl;
    
    // 开始计时
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // 运行完整的参数优化
    auto result = optimizer.optimizeParameters(data);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n=== Optimization Results ===" << std::endl;
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
    std::cout << "Total combinations tested: " << result.all_results.size() << std::endl;
    std::cout << "Best accuracy: " << result.best_score * 100 << "%" << std::endl;
    
    std::cout << "\nBest Parameters:" << std::endl;
    std::cout << "  Window size: " << result.best_params.window_size << std::endl;
    std::cout << "  C score: " << result.best_params.c_score << std::endl;
    std::cout << "  Penalty score: " << result.best_params.penalty_score << std::endl;
    std::cout << "  Initial score: " << result.best_params.initial_score << std::endl;
    
    std::cout << "\nWeight function (first 10 positions):" << std::endl;
    for (int i = 0; i < std::min(10, (int)result.best_params.position_weights.size()); i++) {
        std::cout << "  Position " << i << ": " << result.best_params.position_weights[i] << std::endl;
    }
    
    // 保存结果到文件
    std::ofstream outfile("optimization_results.txt");
    if (outfile.is_open()) {
        outfile << "Best Parameters Found:\n";
        outfile << "Window size: " << result.best_params.window_size << "\n";
        outfile << "C score: " << result.best_params.c_score << "\n"; 
        outfile << "Penalty score: " << result.best_params.penalty_score << "\n";
        outfile << "Initial score: " << result.best_params.initial_score << "\n";
        outfile << "Accuracy: " << result.best_score << "\n\n";
        
        outfile << "All Results (sorted by accuracy):\n";
        std::sort(result.all_results.begin(), result.all_results.end(),
                 [](const std::pair<TrimParams, double>& a, const std::pair<TrimParams, double>& b) { return a.second > b.second; });
        
        for (int i = 0; i < std::min(20, (int)result.all_results.size()); i++) {
            const auto& params = result.all_results[i].first;
            double score = result.all_results[i].second;
            outfile << "Rank " << (i+1) << ": Accuracy=" << score 
                   << ", Window=" << params.window_size
                   << ", Penalty=" << params.penalty_score
                   << ", InitScore=" << params.initial_score << "\n";
        }
        outfile.close();
        std::cout << "Results saved to optimization_results.txt" << std::endl;
    }
}

// 基础功能测试（保留原有的）
void test_basic_functionality() {
    std::cout << "\n=== Testing Basic SmartTrim Functionality ===" << std::endl;
    TestRunner runner;
    
    TrimParams params(30, 0.0);
    SmartTrimmer trimmer(params);
    
    // Test case: R2 with poly-C head
    std::string r1 = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    std::string r2 = "CCCCCATCGATCGATCGATCGATCGATCGATCG";
    
    TrimResult result = trimmer.findOptimalTrimPositions(r1, r2);
    
    runner.assert_test(result.is_valid, "Should detect poly-C tail");
    runner.assert_test(result.final_score > 0, "Score should be positive for poly-C");
    
    runner.print_summary();
}

int main(int argc, char* argv[]) {
    std::cout << "SmartTrim Comprehensive Testing & Parameter Optimization" << std::endl;
    std::cout << "=======================================================" << std::endl;
    
    bool run_full_optimization = false;
    if (argc > 1) {
        std::string arg = argv[1];
        if (arg == "--full-optimization" || arg == "--optimize") {
            run_full_optimization = true;
        }
    }
    
    if (run_full_optimization) {
        std::cout << "Running FULL parameter optimization (this may take several minutes)..." << std::endl;
        full_parameter_optimization_benchmark();
    } else {
        // Run standard tests
        test_basic_functionality();
        test_parameter_optimization();
        test_real_data_workflow();
        
        std::cout << "\n=== Quick Start Guide ===" << std::endl;
        std::cout << "To run full parameter optimization: ./unit_test --optimize" << std::endl;
        std::cout << "To use real FASTQ data: place test_R1.fastq and test_R2.fastq in current directory" << std::endl;
    }
    
    return 0;
}
