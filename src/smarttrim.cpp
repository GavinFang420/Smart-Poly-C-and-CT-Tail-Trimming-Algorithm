#include "smarttrim.h"
#include <algorithm>
#include <iostream>
#include <cmath>

double SmartTrimmer::calculateProgressiveScore(const std::string& merged_seq) {
    double score = params.initial_score;
    int seq_len = merged_seq.length();
    int window_len = std::min(params.window_size, seq_len);
    
    // Focus on the last window_size bases (tail region)
    int start_pos = seq_len - window_len;
    
    for (int i = 0; i < window_len; i++) {
        char base = merged_seq[start_pos + i];
        double position_weight = params.position_weights[i];
        
        if (base == 'C' || base == 'c') {
            score += params.c_score * position_weight;
        } else {
            // A, G, T get penalty (including T, since we focus on polyC not polyC/T)
            score += params.penalty_score * position_weight;
        }
    }
    
    return score;
}

std::pair<int, int> SmartTrimmer::mapToOriginalPositions(int merge_pos, int r1_len, int r2_len) {
    // merge_pos is the position from the end of merged sequence to cut
    // We need to map this back to original R1 and R2 positions
    
    int r1_trim = 0;  // How many bases to trim from R1 tail
    int r2_trim = 0;  // How many bases to trim from R2 head
    
    if (merge_pos > r2_len) {
        // Cut extends into R1
        r2_trim = r2_len;  // Trim entire R2
        r1_trim = merge_pos - r2_len;  // Remaining comes from R1
    } else {
        // Cut only affects R2
        r2_trim = merge_pos;
        r1_trim = 0;
    }
    
    return std::make_pair(r1_trim, r2_trim);
}

TrimResult SmartTrimmer::findOptimalTrimPositions(
    const std::string& r1_seq, 
    const std::string& r2_seq
) {
    TrimResult result;
    
    if (r1_seq.empty() || r2_seq.empty()) {
        return result;  // Invalid input
    }
    
    // Step 1: Merge R1 tail + R2 head (reverse R1 to get tail first)
    std::string r1_tail = r1_seq.substr(std::max(0, (int)r1_seq.length() - params.window_size));
    std::string r2_head = r2_seq.substr(0, std::min(params.window_size, (int)r2_seq.length()));
    
    // Reverse R1 tail so we have: R1_tail_reversed + R2_head
    std::reverse(r1_tail.begin(), r1_tail.end());
    std::string merged = r1_tail + r2_head;
    
    double best_score = params.initial_score;
    int best_cut_pos = 0;
    
    // Step 2: Try different cut positions and find the one with highest score
    for (int cut_pos = 1; cut_pos <= merged.length() && cut_pos <= params.window_size; cut_pos++) {
        // Get the tail segment that would be removed
        std::string tail_segment = merged.substr(merged.length() - cut_pos);
        double score = calculateProgressiveScore(tail_segment);
        
        if (score > best_score) {
            best_score = score;
            best_cut_pos = cut_pos;
        }
    }
    
    // Step 3: Apply cutoff threshold (fixed at 0)
    if (best_score > 0.0) {
        // Map back to original R1/R2 positions
        auto positions = mapToOriginalPositions(best_cut_pos, r1_tail.length(), r2_head.length());
        result.r1_trim_pos = positions.first;
        result.r2_trim_pos = positions.second;
        result.final_score = best_score;
        result.is_valid = true;
    }
    
    return result;
}

std::pair<std::string, std::string> SmartTrimmer::trimReads(
    const std::string& r1_seq,
    const std::string& r2_seq,
    const TrimResult& result
) {
    if (!result.is_valid) {
        return std::make_pair(r1_seq, r2_seq);  // No trimming
    }
    
    // Trim R1 from the end
    std::string trimmed_r1 = r1_seq;
    if (result.r1_trim_pos > 0 && result.r1_trim_pos < r1_seq.length()) {
        trimmed_r1 = r1_seq.substr(0, r1_seq.length() - result.r1_trim_pos);
    }
    
    // Trim R2 from the start
    std::string trimmed_r2 = r2_seq;
    if (result.r2_trim_pos > 0 && result.r2_trim_pos < r2_seq.length()) {
        trimmed_r2 = r2_seq.substr(result.r2_trim_pos);
    }
    
    return std::make_pair(trimmed_r1, trimmed_r2);
}

std::vector<TrimParams> SmartTrimmer::generateParameterMatrix() {
    std::vector<TrimParams> param_list;
    
    // Test different window sizes
    std::vector<int> window_sizes = {20, 23, 25, 30};
    
    // Test different penalty scores
    std::vector<double> penalty_scores = {-2.0, -3.0, -4.0, -5.0};
    
    // Test different initial scores
    std::vector<double> initial_scores = {0.0, -10.0, -20.0};
    
    for (int ws : window_sizes) {
        for (double penalty : penalty_scores) {
            for (double init_score : initial_scores) {
                TrimParams params(ws, init_score);
                params.penalty_score = penalty;
                param_list.push_back(params);
            }
        }
    }
    
    // Generate different weight decay functions
    for (int ws : window_sizes) {
        // Exponential decay
        TrimParams exp_params(ws, 0.0);
        for (int i = 0; i < ws; i++) {
            exp_params.position_weights[i] = std::exp(2.0 * i / (ws - 1));
        }
        param_list.push_back(exp_params);
        
        // Quadratic decay
        TrimParams quad_params(ws, 0.0);
        for (int i = 0; i < ws; i++) {
            double x = (double)i / (ws - 1);
            quad_params.position_weights[i] = 1.0 + 2.0 * x * x;
        }
        param_list.push_back(quad_params);
    }
    
    std::cout << "Generated " << param_list.size() << " parameter configurations for testing." << std::endl;
    return param_list;
}
