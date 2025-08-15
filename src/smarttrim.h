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
#include <utility>
#include "read.h"
#include "sequence.h"

using namespace std;

/**
 * Smart Poly-C/CT Tail Trimming Algorithm
 * 
 * This algorithm addresses the limitations of traditional tail trimming methods
 * by using intelligent scoring based on bisulfite conversion properties:
 * - C bases are rare (~0.4%) in converted sequences
 * - Poly-C/CT tails have high C concentration (85%C + 15%T)
 * - Distance decay: bases closer to sequence ends have higher weights
 */
class SmartTrim {
public:
    // Trimming modes
    enum TrimMode {
        POLY_C_MODE = 0,    // Pure poly-C tail trimming (default)
        POLY_CT_MODE = 1    // Poly-C/CT tail trimming (future development)
    };
    
    // Configuration parameters
    struct Config {
        TrimMode mode = POLY_C_MODE;    // Trimming mode: C-only or C+T
        int window_size = 30;           // Analysis window size (last 30bp)
        int c_score = 10;              // Score for C base appearance
        int t_score = -5;              // Score for T base (varies by mode)
        int penalty_score = -20;       // Penalty for A/G bases
        double end_weight_multiplier = 3.0;  // Weight multiplier for last 10bp
        
        // Poly-C mode validation parameters
        double min_c_ratio = 0.80;     // Minimum C ratio for poly-C validation (80%)
        
        // Poly-CT mode validation parameters (future use)
        double min_ct_ratio = 0.85;    // Minimum C+T ratio for validation (85%)
        double max_t_ratio = 0.15;     // Maximum T ratio in C+T content (15%)
        
        bool enable_validation = true;       // Enable composition validation
        bool enable_position_decay = true;   // Enable distance decay weighting
        int min_trim_length = 1;       // Minimum bases to trim
        int max_trim_length = 25;      // Maximum bases to trim
    };

private:
    Config config;
    
    // Internal scoring and analysis methods
    double calculateBaseScore(char base, int position, int total_length);
    double calculatePositionWeight(int position, int total_length);
    vector<double> calculateProgressiveScores(const string& merged_sequence);
    int findOptimalCutPoint(const vector<double>& scores);
    bool validateCTRatio(const string& sequence, int start_pos, int length);
    
    // Sequence manipulation methods
    string mergeReads(const Read* r1, const Read* r2);
    pair<int, int> mapCutPositions(int cut_point, int r1_length, int r2_length);
    
    // Statistics and validation
    struct TrimStats {
        int total_reads_processed = 0;
        int reads_trimmed = 0;
        int avg_trim_length = 0;
        double avg_c_content = 0.0;
        double avg_ct_ratio = 0.0;
        bool validation_passed = false;
    };
    
    TrimStats stats;

public:
    // Constructor
    SmartTrim();
    SmartTrim(const Config& cfg);
    
    // Main trimming interface
    pair<Read*, Read*> trimReads(const Read* r1, const Read* r2);
    
    // Batch processing
    vector<pair<Read*, Read*>> trimReadPairs(const vector<pair<Read*, Read*>>& read_pairs);
    
    // Configuration management
    void setConfig(const Config& cfg);
    Config getConfig() const;
    
    // Mode management
    void setTrimMode(TrimMode mode);
    TrimMode getTrimMode() const;
    
    // Parameter tuning methods
    void setWindowSize(int size);
    void setCScore(int score);
    void setTScore(int score);          // Separate T score setting
    void setPenaltyScore(int score);
    void setEndWeightMultiplier(double multiplier);
    
    // Validation settings by mode
    void setPolyCValidation(bool enable, double min_c_ratio = 0.80);
    void setPolyCTValidation(bool enable, double min_ct_ratio = 0.85, double max_t_ratio = 0.15);
    void setTrimLengthRange(int min_length, int max_length);
    
    // Statistics and reporting
    TrimStats getStatistics() const;
    void resetStatistics();
    string generateReport() const;
    
    // Validation and testing
    bool testAlgorithmAccuracy(const vector<pair<Read*, Read*>>& test_data);
    void benchmarkPerformance(const vector<pair<Read*, Read*>>& benchmark_data);
    
    // Debug and analysis methods
    vector<double> analyzeSequence(const string& sequence);  // For debugging scoring
    string getDebugInfo(const Read* r1, const Read* r2);    // Detailed analysis output
    
    // Utility methods
    static double calculateCContent(const string& sequence);
    static double calculateCTRatio(const string& sequence);
    static bool isPolyC(const string& sequence, double threshold = 0.8);
    static bool isPolyCT(const string& sequence, double ct_threshold = 0.85, double t_threshold = 0.15);
    
    // Mode-specific analysis
    bool validateSequenceComposition(const string& sequence, TrimMode mode);
    
    // Version and algorithm info
    static string getVersion();
    static string getAlgorithmDescription();
};

// Helper structures for advanced analysis
struct SequenceAnalysis {
    string sequence;
    vector<double> position_scores;
    vector<double> cumulative_scores;
    int optimal_cut_point;
    double c_content;
    double ct_ratio;
    bool trim_recommended;
    string analysis_notes;
};

// Exception classes for error handling
class SmartTrimException : public exception {
private:
    string message;
public:
    SmartTrimException(const string& msg) : message(msg) {}
    virtual const char* what() const noexcept override {
        return message.c_str();
    }
};

class InvalidConfigException : public SmartTrimException {
public:
    InvalidConfigException(const string& msg) : SmartTrimException("Invalid configuration: " + msg) {}
};

class InvalidSequenceException : public SmartTrimException {
public:
    InvalidSequenceException(const string& msg) : SmartTrimException("Invalid sequence: " + msg) {}
};

#endif // SMARTTRIM_H
