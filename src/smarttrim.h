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

// Forward declaration
struct MergeResult;

struct TrimParams {
    int window_size = 30;
    double c_score = 10.0;         // Fixed C score = +10
    double penalty_score = -3.0;   // Penalty for A, G, T
    double initial_score = 0.0;    // Optional initial score
    std::vector<double> position_weights;
    
    // Constructor to initialize position weights
    TrimParams(int ws = 30, double init_score = 0.0) 
        : window_size(ws), initial_score(init_score) {
        // Generate distance decay weights (higher weight toward end)
        position_weights.resize(window_size);
        for (int i = 0; i < window_size; i++) {
            // Linear decay: weight increases from 1.0 to 3.0
            position_weights[i] = 1.0 + 2.0 * (double)i / (window_size - 1);
        }
    }
};

struct TrimResult {
    int r1_trim_pos;    // Position to trim R1 (from end)
    int r2_trim_pos;    // Position to trim R2 (from start)
    double final_score;
    bool is_valid;
    
    TrimResult() : r1_trim_pos(0), r2_trim_pos(0), final_score(0.0), is_valid(false) {}
};

class SmartTrimmer {
private:
    TrimParams params;
    
    // Legacy functions (keep for backward compatibility)
    double calculateProgressiveScore(const std::string& merged_seq);
    std::pair<int, int> mapToOriginalPositions(int merge_pos, int r1_len, int r2_len);

    // New functions for correct merge-based algorithm
    double calculateProgressiveScore(const std::string& sequence, int start_pos, int end_pos);
    TrimResult findOptimalTrimPositions_Window(const std::string& r1_seq, const std::string& r2_seq);
    std::pair<int, int> mapMergedPositionToOriginal(int cut_length_in_r2, 
                                                   const MergeResult& merge_result,
                                                   int original_r1_length, 
                                                   int original_r2_length);

public:
    SmartTrimmer(const TrimParams& p) : params(p) {}
    
    // Main function to find optimal trim positions
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