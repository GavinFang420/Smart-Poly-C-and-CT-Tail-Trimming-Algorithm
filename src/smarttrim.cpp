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

#include "smarttrim.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;

// Constructor implementations
SmartTrim::SmartTrim() {
    // Use default configuration
}

SmartTrim::SmartTrim(const Config& cfg) : config(cfg) {
    // Validate configuration parameters
    if (config.window_size < 10 || config.window_size > 50) {
        throw InvalidConfigException("Window size must be between 10 and 50");
    }
    if (config.c_score <= 0) {
        throw InvalidConfigException("C score must be positive");
    }
    if (config.penalty_score >= 0) {
        throw InvalidConfigException("Penalty score must be negative");
    }
    if (config.min_ct_ratio < 0.5 || config.min_ct_ratio > 1.0) {
        throw InvalidConfigException("Minimum C+T ratio must be between 0.5 and 1.0");
    }
}

// Core algorithm implementation
pair<Read*, Read*> SmartTrim::trimReads(const Read* r1, const Read* r2) {
    if (!r1 || !r2) {
        throw InvalidSequenceException("Null read pointers provided");
    }
    
    stats.total_reads_processed++;
    
    // Step 1: Merge reads and focus on last window_size bp
    string merged_sequence = mergeReads(r1, r2);
    
    if (merged_sequence.length() < config.window_size) {
        // Sequence too short, return original reads
        return make_pair(new Read(*r1), new Read(*r2));
    }
    
    // Step 2: Calculate progressive scores for the analysis window
    string analysis_window = merged_sequence.substr(merged_sequence.length() - config.window_size);
    vector<double> scores = calculateProgressiveScores(analysis_window);
    
    // Step 3: Find optimal cut point
    int cut_point = findOptimalCutPoint(scores);
    
    if (cut_point == 0) {
        // No trimming recommended
        return make_pair(new Read(*r1), new Read(*r2));
    }
    
    // Step 4: Validate C+T ratio if enabled
    if (config.enable_ct_validation) {
        string trimmed_region = analysis_window.substr(analysis_window.length() - cut_point);
        if (!validateCTRatio(trimmed_region, 0, trimmed_region.length())) {
            // Validation failed, return original reads
            return make_pair(new Read(*r1), new Read(*r2));
        }
    }
    
    // Step 5: Map cut positions back to original R1 and R2
    pair<int, int> cut_positions = mapCutPositions(cut_point, r1->length(), r2->length());
    int r1_cut = cut_positions.first;
    int r2_cut = cut_positions.second;
    
    // Step 6: Create trimmed reads
    Read* trimmed_r1 = new Read(*r1);
    Read* trimmed_r2 = new Read(*r2);
    
    // Trim R1 from the end if needed
    if (r1_cut > 0 && r1_cut < r1->length()) {
        string new_seq = r1->mSeq->mStr.substr(0, r1->length() - r1_cut);
        string new_qual = r1->mQuality.substr(0, r1->length() - r1_cut);
        trimmed_r1->mSeq = new Sequence(new_seq);
        trimmed_r1->mQuality = new_qual;
    }
    
    // Trim R2 from the beginning if needed
    if (r2_cut > 0 && r2_cut < r2->length()) {
        string new_seq = r2->mSeq->mStr.substr(r2_cut);
        string new_qual = r2->mQuality.substr(r2_cut);
        trimmed_r2->mSeq = new Sequence(new_seq);
        trimmed_r2->mQuality = new_qual;
    }
    
    // Update statistics
    stats.reads_trimmed++;
    stats.avg_trim_length = (stats.avg_trim_length * (stats.reads_trimmed - 1) + cut_point) / stats.reads_trimmed;
    
    return make_pair(trimmed_r1, trimmed_r2);
}

// Base scoring with position weights
double SmartTrim::calculateBaseScore(char base, int position, int total_length) {
    double base_score = 0.0;
    
    // Base scoring: C gets positive score, others get penalty
    switch (base) {
        case 'C':
        case 'c':
            base_score = config.c_score;
            break;
        case 'T':
        case 't':
            // T gets smaller penalty in CT trimming mode
            base_score = config.penalty_score * 0.3;  // Reduced penalty for T
            break;
        case 'A':
        case 'a':
        case 'G':
        case 'g':
            base_score = config.penalty_score;
            break;
        case 'N':
        case 'n':
            base_score = config.penalty_score * 0.5;  // Moderate penalty for N
            break;
        default:
            base_score = config.penalty_score;
            break;
    }
    
    // Apply position weight if enabled
    if (config.enable_position_decay) {
        double position_weight = calculatePositionWeight(position, total_length);
        base_score *= position_weight;
    }
    
    return base_score;
}

// Position weight calculation with distance decay
double SmartTrim::calculatePositionWeight(int position, int total_length) {
    // Position 0 is the start of analysis window (closer to R2 head)
    // Higher positions are closer to sequence ends (R1 tail)
    
    double relative_position = (double)position / total_length;
    
    // Last 10bp (or 1/3 of window) get higher weight
    int high_weight_region = min(10, total_length / 3);
    
    if (position >= total_length - high_weight_region) {
        // End region: higher weight
        return config.end_weight_multiplier;
    } else {
        // Distance decay: linear decay from middle to start
        // Closer to R2 head gets lower weight
        return 1.0 + relative_position;  // Weight ranges from 1.0 to 2.0
    }
}

// Calculate progressive cumulative scores
vector<double> SmartTrim::calculateProgressiveScores(const string& sequence) {
    vector<double> scores(sequence.length(), 0.0);
    double cumulative_score = 0.0;
    
    // Calculate scores from end to beginning (reverse order for tail trimming)
    for (int i = sequence.length() - 1; i >= 0; i--) {
        double base_score = calculateBaseScore(sequence[i], i, sequence.length());
        cumulative_score += base_score;
        scores[i] = cumulative_score;
    }
    
    return scores;
}

// Find optimal cut point based on scores
int SmartTrim::findOptimalCutPoint(const vector<double>& scores) {
    if (scores.empty()) return 0;
    
    int best_cut_point = 0;
    double best_score = scores[0];  // No trimming score
    
    // Find the position with highest cumulative score
    for (int i = 1; i < scores.size(); i++) {
        if (scores[i] > best_score && i >= config.min_trim_length && i <= config.max_trim_length) {
            best_score = scores[i];
            best_cut_point = scores.size() - i;  // Convert to trim length
        }
    }
    
    // Only trim if score is positive (beneficial)
    return (best_score > 0) ? best_cut_point : 0;
}

// Validate C+T ratio in trimmed region
bool SmartTrim::validateCTRatio(const string& sequence, int start_pos, int length) {
    if (!config.enable_ct_validation) return true;
    
    string region = sequence.substr(start_pos, length);
    if (region.empty()) return false;
    
    int c_count = 0, t_count = 0, total_count = 0;
    
    for (char base : region) {
        switch (base) {
            case 'C': case 'c': c_count++; total_count++; break;
            case 'T': case 't': t_count++; total_count++; break;
            case 'A': case 'a': case 'G': case 'g': total_count++; break;
            // Skip N bases from validation
        }
    }
    
    if (total_count == 0) return false;
    
    double ct_ratio = (double)(c_count + t_count) / total_count;
    double t_in_ct_ratio = (c_count + t_count > 0) ? (double)t_count / (c_count + t_count) : 0.0;
    
    // Update statistics
    stats.avg_ct_ratio = (stats.avg_ct_ratio * (stats.reads_trimmed) + ct_ratio) / (stats.reads_trimmed + 1);
    stats.avg_c_content = (stats.avg_c_content * (stats.reads_trimmed) + (double)c_count / total_count) / (stats.reads_trimmed + 1);
    
    // Validation: C+T ratio >= min_ct_ratio AND T ratio in C+T <= max_t_ratio
    bool validation_result = (ct_ratio >= config.min_ct_ratio) && (t_in_ct_ratio <= config.max_t_ratio);
    stats.validation_passed = validation_result;
    
    return validation_result;
}

// Merge reads for analysis
string SmartTrim::mergeReads(const Read* r1, const Read* r2) {
    // Simple concatenation: R2 + reverse_complement(R1)
    // Focus on the tail region where poly-C/CT sequences are expected
    
    string r2_seq = r2->mSeq->mStr;
    string r1_seq = r1->mSeq->mStr;
    
    // For tail trimming, we're interested in R2 head + R1 tail
    // Merge as: R2_sequence + R1_sequence (R1 is already in correct orientation)
    return r2_seq + r1_seq;
}

// Map cut positions back to original reads
pair<int, int> SmartTrim::mapCutPositions(int cut_point, int r1_length, int r2_length) {
    // cut_point represents bases to trim from the merged sequence tail
    // We need to determine how many bases to trim from R1 tail and R2 head
    
    int r2_cut = min(cut_point, r2_length);  // Trim from R2 beginning
    int r1_cut = cut_point - r2_cut;         // Remaining trim from R1 end
    r1_cut = max(0, min(r1_cut, r1_length)); // Ensure valid range
    
    return make_pair(r1_cut, r2_cut);
}

// Batch processing
vector<pair<Read*, Read*>> SmartTrim::trimReadPairs(const vector<pair<Read*, Read*>>& read_pairs) {
    vector<pair<Read*, Read*>> trimmed_pairs;
    trimmed_pairs.reserve(read_pairs.size());
    
    for (const auto& read_pair : read_pairs) {
        trimmed_pairs.push_back(trimReads(read_pair.first, read_pair.second));
    }
    
    return trimmed_pairs;
}

// Configuration management
void SmartTrim::setConfig(const Config& cfg) {
    config = cfg;
}

SmartTrim::Config SmartTrim::getConfig() const {
    return config;
}

// Parameter tuning methods
void SmartTrim::setWindowSize(int size) {
    if (size < 10 || size > 50) {
        throw InvalidConfigException("Window size must be between 10 and 50");
    }
    config.window_size = size;
}

void SmartTrim::setCScore(int score) {
    if (score <= 0) {
        throw InvalidConfigException("C score must be positive");
    }
    config.c_score = score;
}

void SmartTrim::setPenaltyScore(int score) {
    if (score >= 0) {
        throw InvalidConfigException("Penalty score must be negative");
    }
    config.penalty_score = score;
}

void SmartTrim::setEndWeightMultiplier(double multiplier) {
    if (multiplier < 1.0) {
        throw InvalidConfigException("End weight multiplier must be >= 1.0");
    }
    config.end_weight_multiplier = multiplier;
}

void SmartTrim::setCTValidation(bool enable, double min_ct_ratio, double max_t_ratio) {
    config.enable_ct_validation = enable;
    config.min_ct_ratio = min_ct_ratio;
    config.max_t_ratio = max_t_ratio;
}

void SmartTrim::setTrimLengthRange(int min_length, int max_length) {
    if (min_length < 0 || max_length < min_length) {
        throw InvalidConfigException("Invalid trim length range");
    }
    config.min_trim_length = min_length;
    config.max_trim_length = max_length;
}

// Statistics and reporting
SmartTrim::TrimStats SmartTrim::getStatistics() const {
    return stats;
}

void SmartTrim::resetStatistics() {
    stats = TrimStats();
}

string SmartTrim::generateReport() const {
    ostringstream report;
    
    report << "=== Smart Poly-C/CT Tail Trimming Report ===" << endl;
    report << "Total reads processed: " << stats.total_reads_processed << endl;
    report << "Reads trimmed: " << stats.reads_trimmed << endl;
    report << "Trimming rate: " << fixed << setprecision(2) 
           << (stats.total_reads_processed > 0 ? (double)stats.reads_trimmed / stats.total_reads_processed * 100 : 0) << "%" << endl;
    report << "Average trim length: " << stats.avg_trim_length << " bp" << endl;
    report << "Average C content: " << fixed << setprecision(4) << stats.avg_c_content * 100 << "%" << endl;
    report << "Average C+T ratio: " << fixed << setprecision(4) << stats.avg_ct_ratio * 100 << "%" << endl;
    report << "Last validation result: " << (stats.validation_passed ? "PASSED" : "FAILED") << endl;
    
    report << "\n=== Algorithm Configuration ===" << endl;
    report << "Window size: " << config.window_size << " bp" << endl;
    report << "C score: " << config.c_score << endl;
    report << "Penalty score: " << config.penalty_score << endl;
    report << "End weight multiplier: " << config.end_weight_multiplier << "x" << endl;
    report << "C+T validation: " << (config.enable_ct_validation ? "ENABLED" : "DISABLED") << endl;
    if (config.enable_ct_validation) {
        report << "  Min C+T ratio: " << fixed << setprecision(2) << config.min_ct_ratio * 100 << "%" << endl;
        report << "  Max T ratio in C+T: " << fixed << setprecision(2) << config.max_t_ratio * 100 << "%" << endl;
    }
    
    return report.str();
}

// Utility methods
double SmartTrim::calculateCContent(const string& sequence) {
    if (sequence.empty()) return 0.0;
    
    int c_count = 0;
    for (char base : sequence) {
        if (base == 'C' || base == 'c') {
            c_count++;
        }
    }
    
    return (double)c_count / sequence.length();
}

double SmartTrim::calculateCTRatio(const string& sequence) {
    if (sequence.empty()) return 0.0;
    
    int ct_count = 0, total_count = 0;
    for (char base : sequence) {
        if (base == 'C' || base == 'c' || base == 'T' || base == 't') {
            ct_count++;
        }
        if (base != 'N' && base != 'n') {
            total_count++;
        }
    }
    
    return total_count > 0 ? (double)ct_count / total_count : 0.0;
}

bool SmartTrim::isPolyC(const string& sequence, double threshold) {
    return calculateCContent(sequence) >= threshold;
}

bool SmartTrim::isPolyCT(const string& sequence, double ct_threshold, double t_threshold) {
    double ct_ratio = calculateCTRatio(sequence);
    if (ct_ratio < ct_threshold) return false;
    
    // Check T ratio within C+T content
    int c_count = 0, t_count = 0;
    for (char base : sequence) {
        if (base == 'C' || base == 'c') c_count++;
        else if (base == 'T' || base == 't') t_count++;
    }
    
    if (c_count + t_count == 0) return false;
    
    double t_in_ct_ratio = (double)t_count / (c_count + t_count);
    return t_in_ct_ratio <= t_threshold;
}

// Debug and analysis methods
vector<double> SmartTrim::analyzeSequence(const string& sequence) {
    return calculateProgressiveScores(sequence);
}

string SmartTrim::getDebugInfo(const Read* r1, const Read* r2) {
    ostringstream debug;
    
    string merged = mergeReads(r1, r2);
    string analysis_window = merged.substr(max(0, (int)merged.length() - config.window_size));
    vector<double> scores = calculateProgressiveScores(analysis_window);
    int cut_point = findOptimalCutPoint(scores);
    
    debug << "=== SmartTrim Debug Info ===" << endl;
    debug << "R1 length: " << r1->length() << " bp" << endl;
    debug << "R2 length: " << r2->length() << " bp" << endl;
    debug << "Merged length: " << merged.length() << " bp" << endl;
    debug << "Analysis window: " << analysis_window << endl;
    debug << "Recommended cut point: " << cut_point << " bp" << endl;
    debug << "C content in window: " << fixed << setprecision(2) 
          << calculateCContent(analysis_window) * 100 << "%" << endl;
    debug << "C+T ratio in window: " << fixed << setprecision(2) 
          << calculateCTRatio(analysis_window) * 100 << "%" << endl;
    
    return debug.str();
}

// Version and algorithm info
string SmartTrim::getVersion() {
    return "SmartTrim v1.0.0 - Smart Poly-C/CT Tail Trimming Algorithm";
}

string SmartTrim::getAlgorithmDescription() {
    return "Smart trimming algorithm based on bisulfite conversion properties. "
           "Uses intelligent scoring with distance decay to identify optimal poly-C/CT regions. "
           "Preserves valuable R1 data while accurately trimming artificial tails.";
}
