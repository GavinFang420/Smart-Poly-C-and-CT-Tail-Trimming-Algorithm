#include "smarttrim.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <chrono>

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

void test_basic_functionality() {
    std::cout << "\n=== Testing Basic Functionality ===" << std::endl;
    TestRunner runner;
    
    TrimParams params(30, 0.0);
    SmartTrimmer trimmer(params);
    
    // Test case 1: R2 with poly-C head
    std::string r1 = "ATCGATCGATCGATCGATCGATCGATCGATCG";  // 32bp, normal sequence
    std::string r2 = "CCCCCATCGATCGATCGATCGATCGATCGATCG";  // 32bp, 5 C's at start
    
    TrimResult result = trimmer.findOptimalTrimPositions(r1, r2);
    
    runner.assert_test(result.is_valid, "Should detect poly-C tail");
    runner.assert_test(result.final_score > 0, "Score should be positive for poly-C");
    runner.assert_test(result.r2_trim_pos > 0, "Should trim R2");
    
    // Test case 2: No poly-C, should not trim
    r1 = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    r2 = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    
    result = trimmer.findOptimalTrimPositions(r1, r2);
    
    runner.assert_test(!result.is_valid || result.final_score <= 0, "Should not trim normal sequences");
    
    // Test case 3: R1 with C-rich tail
    r1 = "ATCGATCGATCGATCGATCGATCGATCCCCC";  // C's at R1 tail
    r2 = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    
    result = trimmer.findOptimalTrimPositions(r1, r2);
    
    runner.assert_test(result.is_valid, "Should detect C-rich R1 tail");
    runner.assert_test(result.r1_trim_pos > 0, "Should trim R1 tail");
    
    runner.print_summary();
}

void test_trimming_functionality() {
    std::cout << "\n=== Testing Trimming Functionality ===" << std::endl;
    TestRunner runner;
    
    TrimParams params(30, 0.0);
    SmartTrimmer trimmer(params);
    
    std::string r1 = "ATCGATCGATCGATCGATCGATCGATCCCCC";  // 30bp
    std::string r2 = "CCCCCATCGATCGATCGATCGATCGATCGATCG";  // 32bp
    
    TrimResult result = trimmer.findOptimalTrimPositions(r1, r2);
    
    if (result.is_valid) {
        auto trimmed = trimmer.trimReads(r1, r2, result);
        
        runner.assert_test(trimmed.first.length() < r1.length() || 
                          trimmed.second.length() < r2.length(), 
                          "At least one read should be trimmed");
        
        runner.assert_test(trimmed.first.length() > 0 && trimmed.second.length() > 0, 
                          "Trimmed reads should not be empty");
        
        std::cout << "Original R1 length: " << r1.length() << " -> " << trimmed.first.length() << std::endl;
        std::cout << "Original R2 length: " << r2.length() << " -> " << trimmed.second.length() << std::endl;
    }
    
    runner.print_summary();
}

void test_parameter_matrix() {
    std::cout << "\n=== Testing Parameter Matrix Generation ===" << std::endl;
    TestRunner runner;
    
    auto param_matrix = SmartTrimmer::generateParameterMatrix();
    
    runner.assert_test(param_matrix.size() > 50, "Should generate substantial number of parameter configs");
    runner.assert_test(param_matrix.size() < 200, "Should not generate excessive parameters");
    
    // Verify different configurations
    bool has_different_window_sizes = false;
    bool has_different_penalties = false;
    
    for (size_t i = 1; i < param_matrix.size(); i++) {
        if (param_matrix[i].window_size != param_matrix[0].window_size) {
            has_different_window_sizes = true;
        }
        if (param_matrix[i].penalty_score != param_matrix[0].penalty_score) {
            has_different_penalties = true;
        }
    }
    
    runner.assert_test(has_different_window_sizes, "Should have different window sizes");
    runner.assert_test(has_different_penalties, "Should have different penalty scores");
    
    runner.print_summary();
}

void test_edge_cases() {
    std::cout << "\n=== Testing Edge Cases ===" << std::endl;
    TestRunner runner;
    
    TrimParams params(30, 0.0);
    SmartTrimmer trimmer(params);
    
    // Test empty strings
    TrimResult result = trimmer.findOptimalTrimPositions("", "ATCG");
    runner.assert_test(!result.is_valid, "Should handle empty R1");
    
    result = trimmer.findOptimalTrimPositions("ATCG", "");
    runner.assert_test(!result.is_valid, "Should handle empty R2");
    
    // Test very short sequences
    result = trimmer.findOptimalTrimPositions("AT", "CG");
    runner.assert_test(true, "Should handle very short sequences without crashing");
    
    // Test sequences shorter than window size
    result = trimmer.findOptimalTrimPositions("ATCGATCG", "CCCCCCC");
    runner.assert_test(true, "Should handle sequences shorter than window size");
    
    // Test all C sequence (extreme case)
    result = trimmer.findOptimalTrimPositions("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", 
                                             "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    runner.assert_test(result.is_valid, "Should handle all-C sequences");
    runner.assert_test(result.final_score > 0, "All-C should get positive score");
    
    runner.print_summary();
}

void benchmark_performance() {
    std::cout << "\n=== Performance Benchmark ===" << std::endl;
    
    TrimParams params(30, 0.0);
    SmartTrimmer trimmer(params);
    
    // Generate test data
    std::string r1(150, 'A');  // 150bp read
    std::string r2(150, 'A');
    
    // Add some C's to make it interesting
    for (int i = 140; i < 150; i++) {
        r1[i] = 'C';
    }
    for (int i = 0; i < 10; i++) {
        r2[i] = 'C';
    }
    
    const int num_tests = 10000;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_tests; i++) {
        TrimResult result = trimmer.findOptimalTrimPositions(r1, r2);
        if (result.is_valid) {
            trimmer.trimReads(r1, r2, result);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    double avg_time = (double)duration.count() / num_tests;
    
    std::cout << "Processed " << num_tests << " read pairs in " 
              << duration.count() << " microseconds" << std::endl;
    std::cout << "Average time per read pair: " << avg_time << " microseconds" << std::endl;
    std::cout << "Estimated throughput: " << (int)(1000000.0 / avg_time) 
              << " read pairs per second" << std::endl;
    
    if (avg_time < 100.0) {
        std::cout << "[PASS] Performance is acceptable (< 100μs per read pair)" << std::endl;
    } else {
        std::cout << "[WARN] Performance might need optimization (> 100μs per read pair)" << std::endl;
    }
}

void run_parameter_matrix_test() {
    std::cout << "\n=== Running Parameter Matrix Test (100+ configurations) ===" << std::endl;
    
    auto param_matrix = SmartTrimmer::generateParameterMatrix();
    
    // Test data with clear poly-C pattern
    std::string r1 = "ATCGATCGATCGATCGATCGATCGATCCCCCCCCCCC";  // C's at end
    std::string r2 = "CCCCCCCCCATCGATCGATCGATCGATCGATCGATCG";   // C's at start
    
    int configs_with_trimming = 0;
    double best_score = -1000.0;
    TrimParams best_params(30, 0.0);
    
    std::cout << "Testing " << param_matrix.size() << " parameter configurations..." << std::endl;
    
    for (size_t i = 0; i < param_matrix.size(); i++) {
        SmartTrimmer trimmer(param_matrix[i]);
        TrimResult result = trimmer.findOptimalTrimPositions(r1, r2);
        
        if (result.is_valid && result.final_score > 0) {
            configs_with_trimming++;
            if (result.final_score > best_score) {
                best_score = result.final_score;
                best_params = param_matrix[i];
            }
        }
        
        if ((i + 1) % 20 == 0) {
            std::cout << "Processed " << (i + 1) << " configurations..." << std::endl;
        }
    }
    
    std::cout << "\n=== Parameter Matrix Results ===" << std::endl;
    std::cout << "Configurations that detected poly-C: " << configs_with_trimming 
              << " / " << param_matrix.size() << std::endl;
    std::cout << "Best score achieved: " << best_score << std::endl;
    std::cout << "Best parameters:" << std::endl;
    std::cout << "  Window size: " << best_params.window_size << std::endl;
    std::cout << "  Initial score: " << best_params.initial_score << std::endl;
    std::cout << "  Penalty score: " << best_params.penalty_score << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "SmartTrim Unit Tests" << std::endl;
    std::cout << "===================" << std::endl;
    
    // Check for benchmark flag
    bool run_benchmark = false;
    if (argc > 1 && std::string(argv[1]) == "--benchmark") {
        run_benchmark = true;
    }
    
    if (run_benchmark) {
        run_parameter_matrix_test();
        benchmark_performance();
    } else {
        // Run standard unit tests
        test_basic_functionality();
        test_trimming_functionality();
        test_parameter_matrix();
        test_edge_cases();
        benchmark_performance();
    }
    
    std::cout << "\nAll tests completed!" << std::endl;
    return 0;
}
